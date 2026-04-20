"""DeconvolutionLab2 (Java) backend, called via PyImageJ.

Mirrors the original MATLAB pipeline verbatim — the same plugin, the same
class, the same arguments::

    result = DL2.RL(image, psf, 30.0, '-out stack short')

This backend is a thin PyImageJ bridge: it starts a single Fiji JVM per
Python process in ``interactive`` mode (AWT enabled, no main window),
resolves the ``DL2`` class from the ``DeconvolutionLab_2.jar`` plugin,
and invokes ``DL2.RL(...)`` on each stack. Image I/O is handled by
ImageJ itself (``IJ.openImage`` / ``IJ.saveAsTiff``), so there is no
numpy <-> Java conversion step.

**Display required.** DL2 internally constructs AWT/Swing components
even for programmatic calls, so the JVM can't be fully headless. This
engine works on a local Mac/Linux desktop but will NOT work in a
server/CI/SSH-only environment without X forwarding. Use the ``scipy``
or ``torch`` engines in those cases.

Because the JVM cannot be restarted inside one process, the gateway is
cached in a module global and reused. Separate ``DL2.RL`` calls on
separate data are safe to run on multiple Java threads, which is what
the CLI's channel-level parallelism relies on.
"""

from __future__ import annotations

import atexit
import logging
import os
import threading
from pathlib import Path

log = logging.getLogger(__name__)


_GATEWAY_LOCK = threading.Lock()
_IJ_GATEWAY = None  # pyimagej gateway (initialized lazily)
_DL2_CLASS = None
_IJ_CLASS = None


class DL2NotAvailableError(RuntimeError):
    """Raised when the DL2 engine can't be loaded (missing Fiji/plugin/JVM)."""


def _resolve_fiji_dir(fiji_dir: Path | str | None) -> Path:
    """Find a Fiji.app directory, preferring an explicit argument, then env var."""
    candidates: list[Path] = []
    if fiji_dir is not None:
        candidates.append(Path(fiji_dir))
    env_dir = os.environ.get("FIJI_DIR")
    if env_dir:
        candidates.append(Path(env_dir))
    for cand in candidates:
        if cand.is_dir():
            return cand.resolve()
    raise DL2NotAvailableError(
        "No Fiji installation found for the DL2 engine.\n"
        "Either pass --fiji-dir /path/to/Fiji.app, set FIJI_DIR in your "
        "environment, or use --engine scipy/torch to skip DL2."
    )


def _verify_dl2_plugin(fiji_dir: Path) -> None:
    """Sanity-check that DeconvolutionLab_2.jar is actually on the plugin path."""
    plugins = fiji_dir / "plugins"
    if not plugins.is_dir():
        raise DL2NotAvailableError(
            f"Fiji directory has no plugins/ subfolder: {fiji_dir}"
        )
    matches = list(plugins.glob("DeconvolutionLab*.jar")) + list(
        plugins.glob("**/DeconvolutionLab*.jar")
    )
    if not matches:
        raise DL2NotAvailableError(
            f"DeconvolutionLab_2.jar not found under {plugins}.\n"
            "Install the plugin via Fiji > Help > Update... > Manage Update Sites,\n"
            "enable 'DeconvolutionLab2', then run the updater."
        )


def _get_gateway(fiji_dir: Path | str | None):
    """Lazy-init a single shared PyImageJ gateway + DL2 / IJ class handles."""
    global _IJ_GATEWAY, _DL2_CLASS, _IJ_CLASS
    with _GATEWAY_LOCK:
        if _IJ_GATEWAY is not None:
            return _IJ_GATEWAY, _DL2_CLASS, _IJ_CLASS
        try:
            import imagej
            import scyjava
        except ImportError as err:
            raise DL2NotAvailableError(
                "The DL2 engine requires pyimagej and scyjava.\n"
                "Install with: pip install pyimagej scyjava\n"
                "You also need a JDK (e.g. `brew install openjdk` on macOS)."
            ) from err

        resolved = _resolve_fiji_dir(fiji_dir)
        _verify_dl2_plugin(resolved)

        # DL2 internally constructs AWT/Swing components (progress monitors,
        # etc.) even for programmatic ``DL2.RL`` calls. PyImageJ's
        # ``mode="headless"`` forces ``java.awt.headless=true``, which makes
        # any AWT call throw ``HeadlessException``. ``mode="interactive"``
        # starts the JVM with AWT enabled but doesn't pop the Fiji main
        # window — good enough for DL2, as long as the user has a display
        # (local Mac/Linux desktop). This will NOT work over an SSH session
        # without X forwarding, or in a server/CI environment.
        log.info("Starting Fiji JVM (interactive mode for DL2 AWT) from %s ...", resolved)
        ij = imagej.init(str(resolved), mode="interactive")
        log.info("Fiji %s ready.", ij.getVersion())

        try:
            DL2 = scyjava.jimport("DL2")
        except Exception as err:
            raise DL2NotAvailableError(
                "Fiji started but the DL2 class could not be loaded. "
                "Is DeconvolutionLab_2.jar in the plugins/ folder?"
            ) from err

        IJ = scyjava.jimport("ij.IJ")

        _IJ_GATEWAY = ij
        _DL2_CLASS = DL2
        _IJ_CLASS = IJ

        atexit.register(_dispose_gateway)
        return ij, DL2, IJ


def _dispose_gateway() -> None:
    """Best-effort JVM shutdown on interpreter exit."""
    global _IJ_GATEWAY
    if _IJ_GATEWAY is not None:
        try:
            _IJ_GATEWAY.dispose()
        except Exception:  # pragma: no cover — interpreter-shutdown path
            pass
        _IJ_GATEWAY = None


def deconvolve_file(
    input_path: Path,
    output_path: Path,
    psf_path: Path,
    num_iter: int,
    fiji_dir: Path | str | None = None,
) -> Path:
    """Deconvolve a single multi-page TIFF via ``DL2.RL``.

    Uses ImageJ's own TIFF I/O: reads ``input_path`` + ``psf_path`` with
    ``IJ.openImage`` (so ImageJ handles the stack dimensionality the same
    way the MATLAB pipeline did), runs ``DL2.RL(image, psf, num_iter,
    '-out stack short')``, and saves the result with ``IJ.saveAsTiff``
    (uncompressed, uint16 — matching the MATLAB output).

    The ``-out stack short`` option mirrors the MATLAB argument verbatim.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    psf_path = Path(psf_path)

    _, DL2, IJ = _get_gateway(fiji_dir)

    imp_image = IJ.openImage(str(input_path))
    if imp_image is None:
        raise RuntimeError(f"IJ.openImage returned null for image: {input_path}")

    imp_psf = IJ.openImage(str(psf_path))
    if imp_psf is None:
        imp_image.close()
        raise RuntimeError(f"IJ.openImage returned null for PSF: {psf_path}")

    try:
        # '-out stack short' matches the MATLAB output type (uint16 ImagePlus).
        result_imp = DL2.RL(imp_image, imp_psf, float(num_iter),
                            "-out stack short")
        if result_imp is None:
            raise RuntimeError(
                f"DL2.RL returned null for {input_path.name} "
                "(check Fiji console; common causes: image/PSF dimensionality "
                "mismatch, unreadable PSF)."
            )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        IJ.saveAsTiff(result_imp, str(output_path))
        try:
            result_imp.close()
        except Exception:
            pass
    finally:
        try:
            imp_image.close()
        except Exception:
            pass
        try:
            imp_psf.close()
        except Exception:
            pass

    return output_path
