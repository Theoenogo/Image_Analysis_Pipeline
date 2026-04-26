"""DeconvolutionLab2 backend — runs DL2 via Fiji CLI as a subprocess.

Mirrors the original MATLAB pipeline: same plugin, same algorithm, same
``-out stack short`` output type. Instead of embedding the JVM inside
the Python process (which fails on macOS due to AWT/AppKit constraints),
this backend launches Fiji as a normal macOS/Linux application via
``subprocess``, hands it a tiny ImageJ macro that calls DL2, and picks
up the output TIFF when Fiji exits.

Each ``deconvolve_file`` call spawns one Fiji process. Channel-level
parallelism (GFP + Cy) is handled by the caller launching two
``deconvolve_file`` calls concurrently — they run as separate OS
processes with no shared JVM state.
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from pathlib import Path

log = logging.getLogger(__name__)


class DL2NotAvailableError(RuntimeError):
    """Raised when the DL2 engine can't be used (missing Fiji/plugin)."""


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


def _find_fiji_launcher(fiji_dir: Path) -> Path:
    """Locate the Fiji launcher executable for the current platform."""
    candidates = [
        fiji_dir / "Contents" / "MacOS" / "ImageJ-macosx",
        fiji_dir / "ImageJ-linux64",
        fiji_dir / "ImageJ-win64.exe",
    ]
    for c in candidates:
        if c.is_file():
            return c
    raise DL2NotAvailableError(
        f"Cannot find Fiji launcher in {fiji_dir}.\n"
        f"Tried: {', '.join(c.name for c in candidates)}"
    )


def _verify_dl2_plugin(fiji_dir: Path) -> None:
    """Sanity-check that DeconvolutionLab_2.jar is on the plugin path."""
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


def deconvolve_file(
    input_path: Path,
    output_path: Path,
    psf_path: Path,
    num_iter: int,
    fiji_dir: Path | str | None = None,
) -> Path:
    """Deconvolve a single multi-page TIFF via Fiji CLI + DL2 macro.

    Writes a temporary ImageJ macro that calls DL2's ``Run`` command,
    launches Fiji as a subprocess, and picks up the output TIFF.
    Fiji runs as a normal application (not headless) so DL2's internal
    AWT/Swing components work without issue.
    """
    input_path = Path(input_path).resolve()
    output_path = Path(output_path).resolve()
    psf_path = Path(psf_path).resolve()

    resolved_fiji = _resolve_fiji_dir(fiji_dir)
    _verify_dl2_plugin(resolved_fiji)
    launcher = _find_fiji_launcher(resolved_fiji)

    # ImageJ macros use forward slashes on all platforms. Windows
    # backslashes inside macro strings break path parsing.
    def _macro_path(p: Path) -> str:
        return str(p).replace("\\", "/")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Load image + PSF via ImageJ's own opener, then hand them to DL2
    # by title using '-image platform' / '-psf platform'. DL2's built-in
    # '-image file ...' opener rejects some ImageJ-hyperstack TIFFs with
    # "Image: Not valid", but IJ.openImage() handles them fine.
    # -out stack short NAME writes NAME.tif (uint16) into -path dir.
    # NB: DL2 only strips bracket-quoting from its -image/-psf file args;
    # wrapping -path in [...] makes DL2 treat the brackets as literal
    # characters and the directory is reported as "(not-valid)", which
    # silently skips the output save. So -path must be bare.
    # NB: DL2's "Run" command dispatches the actual RL onto a background
    # worker thread and returns to the macro immediately. If the macro
    # calls run("Quit") right after, Fiji shuts down before the worker
    # finishes and the output file is never written. So after kicking
    # DL2 off we poll for the expected output TIFF (up to poll_seconds)
    # and only quit once it exists. The outer subprocess timeout still
    # bounds total wall-clock time.
    poll_seconds = 580  # a bit under subprocess timeout (600s) so we
                        # can log a clean timeout before Python kills us.
    img_title = "dl2_image"
    psf_title = "dl2_psf"
    expected_output = output_path.parent / f"{output_path.stem}.tif"
    m_input = _macro_path(input_path)
    m_psf = _macro_path(psf_path)
    m_outfile = _macro_path(expected_output)
    m_outdir = _macro_path(output_path.parent)
    macro = (
        f'print("DL2 starting: {input_path.name}");\n'
        f'open("{m_input}");\n'
        f'rename("{img_title}");\n'
        f'open("{m_psf}");\n'
        f'rename("{psf_title}");\n'
        f'outFile = "{m_outfile}";\n'
        f'run("DeconvolutionLab2 Run", '
        f'"-image platform {img_title} '
        f'-psf platform {psf_title} '
        f'-algorithm RL {num_iter} '
        f'-out stack short {output_path.stem} '
        f'-path {m_outdir}");\n'
        f'print("DL2 Run returned to macro, waiting for output file");\n'
        f'waited = 0;\n'
        f'while (!File.exists(outFile) && waited < {poll_seconds}) {{\n'
        f'    wait(1000);\n'
        f'    waited = waited + 1;\n'
        f'}}\n'
        f'if (File.exists(outFile)) {{\n'
        f'    print("DL2 output appeared after " + waited + "s: " + outFile);\n'
        f'}} else {{\n'
        f'    print("DL2 output never appeared, gave up after " + waited + "s");\n'
        f'}}\n'
        f'print("DL2 finished: {input_path.name}");\n'
        f'run("Quit");\n'
    )
    log.debug("DL2 macro:\n%s", macro)

    fd, macro_path = tempfile.mkstemp(suffix=".ijm")
    try:
        os.write(fd, macro.encode())
        os.close(fd)

        log.info("Running DL2 via Fiji CLI: %s", input_path.name)

        # DL2 hardcodes a JFrame in Deconvolution.deconvolve(), so Fiji
        # cannot run with --headless — the plugin always throws
        # HeadlessException. Launch the Fiji launcher binary directly so
        # AWT is available; on macOS a subprocess spawned from a GUI
        # terminal session inherits display access, so a brief Fiji
        # window appears per deconvolution but the macro quits on its own.
        cmd = [str(launcher), "-macro", macro_path]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,
        )

        if result.stdout:
            for line in result.stdout.strip().splitlines():
                log.info("Fiji: %s", line)
        if result.stderr:
            for line in result.stderr.strip().splitlines():
                log.debug("Fiji stderr: %s", line)
    finally:
        try:
            os.unlink(macro_path)
        except OSError:
            pass

    if output_path.is_file():
        return output_path

    # DL2 may have saved without extension match — search for the stem
    alt = output_path.parent / f"{output_path.stem}.tif"
    if alt.is_file() and alt != output_path:
        alt.rename(output_path)
        return output_path

    recent = sorted(
        output_path.parent.glob("*.tif"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )[:10]
    raise RuntimeError(
        f"DL2 did not produce output at {output_path}.\n"
        f"Return code: {result.returncode}\n"
        f"Recent .tif in output dir: {[f.name for f in recent]}\n"
        f"--- stdout ---\n{result.stdout or '(empty)'}\n"
        f"--- stderr ---\n{result.stderr or '(empty)'}"
    )
