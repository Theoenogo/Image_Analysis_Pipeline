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

import numpy as np
import tifffile

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

        # DL2 internally constructs AWT/Swing components even for programmatic
        # DL2.RL calls, so the JVM must NOT be fully headless.
        #
        # On macOS, PyImageJ's mode="interactive" refuses to start when the
        # CoreFoundation/AppKit event loop isn't already running (which it
        # never is for a plain `python` CLI). mode="interactive:force" skips
        # that guard while still initialising the JVM with AWT enabled, which
        # is what DL2 needs. We pass `-system` to DL2.RL so it writes its
        # progress to stdout instead of trying to open a Swing monitor window,
        # meaning nothing tries to actually render to screen and the lack of
        # an event loop never matters.
        log.info("Starting Fiji JVM (interactive:force for DL2 AWT) from %s ...", resolved)
        ij = imagej.init(str(resolved), mode="interactive:force")
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


def _numpy_to_imageplus(arr: np.ndarray, ij_gateway, title: str = "image"):
    """Convert a (Z, Y, X) numpy array to an ImageJ1 ImagePlus.

    Uses PyImageJ's ``ij.py.to_imageplus()`` conversion, which constructs
    the ImagePlus directly from array data without any file I/O or
    display interaction. This works reliably in the ``interactive:force``
    macOS CLI environment where ``IJ.openImage`` returns pixel-empty
    ImagePlus objects.
    """
    arr = np.asarray(arr)
    if arr.ndim == 2:
        arr = arr[np.newaxis, ...]
    imp = ij_gateway.py.to_imageplus(arr.astype(np.float32))
    if imp is None:
        raise RuntimeError(f"ij.py.to_imageplus returned null for: {title}")
    log.debug("Built ImagePlus '%s': %dx%dx%d slices",
              title, imp.getWidth(), imp.getHeight(), imp.getStackSize())
    return imp


def _load_as_imageplus(path: Path, ij_gateway):
    """Read a TIFF with tifffile and convert to ImagePlus via PyImageJ."""
    arr = tifffile.imread(str(path))
    return _numpy_to_imageplus(arr, ij_gateway, title=path.name)


def deconvolve_file(
    input_path: Path,
    output_path: Path,
    psf_path: Path,
    num_iter: int,
    fiji_dir: Path | str | None = None,
) -> Path:
    """Deconvolve a single multi-page TIFF via ``DL2.RL``.

    Reads ``input_path`` + ``psf_path`` using tifffile + PyImageJ's
    ``to_imageplus()`` — bypassing ``IJ.openImage`` entirely to avoid
    the macOS ``interactive:force`` issue where that call returns a
    pixel-empty ImagePlus. Runs ``DL2.RL`` and saves with
    ``IJ.saveAsTiff`` (uncompressed, uint16 — matching MATLAB output).
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    psf_path = Path(psf_path)

    ij_gateway, DL2, IJ = _get_gateway(fiji_dir)

    imp_image = _load_as_imageplus(input_path, ij_gateway)
    imp_psf = _load_as_imageplus(psf_path, ij_gateway)

    try:
        # '-out stack short' matches the MATLAB output type (uint16 ImagePlus).
        # '-system' redirects DL2's progress log to stdout instead of opening
        # a Swing monitor window — required when no event loop is running.
        result_imp = DL2.RL(imp_image, imp_psf, float(num_iter),
                            "-out stack short -system")
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
