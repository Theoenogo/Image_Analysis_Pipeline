"""CLI entry point: deconvolve a single folder of images for one channel.

Standalone version of the deconvolution step in ``preprocess.py``, for when
you just want to re-run Richardson-Lucy on one folder of TIFF stacks against
one PSF, without the full stacking/consolidation pipeline.

    python deconvolve_folder.py \\
        --input-dir  /path/to/sampleA/gfp \\
        --psf        /path/to/gfp_psf.tif \\
        --output-dir /path/to/sampleA/gfp_decon

Every ``*.tif`` directly inside ``--input-dir`` (non-recursive) is
deconvolved and written to ``<output-dir>/<name>_decon.tif``. If
``--output-dir`` is omitted, outputs go to ``<input-dir>_decon`` next to the
input folder.

Pass ``--jobs N`` to deconvolve multiple images concurrently (default: 1,
sequential). For ``--engine dl2`` each job is its own Fiji process, so this
can meaningfully speed up large batches.

This script reuses the same engines, PSF auto-crop, and output rescaling as
``preprocess.py`` (see ``src/preprocess/deconvolve.py``); it does not modify
any of that code.
"""

from __future__ import annotations

import argparse
import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# Make the src/ layout importable without requiring `pip install -e .`.
_SRC = Path(__file__).parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from preprocess.deconvolve import (  # noqa: E402
    DeconvolutionConfig,
    _autocrop_psf,
    _deconvolve_file_python,
    _load_tiff_stack,
    _save_cropped_psf,
)

log = logging.getLogger("deconvolve_folder")


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="deconvolve_folder",
        description="Run Richardson-Lucy deconvolution on every TIFF stack "
        "in a single folder, against a single PSF.",
    )
    p.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Folder containing .tif stacks to deconvolve (non-recursive).",
    )
    p.add_argument(
        "--psf",
        type=Path,
        required=True,
        help="Measured PSF TIFF for this channel.",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Folder to write <name>_decon.tif outputs to. "
        "Default: <input-dir>_decon next to the input folder.",
    )
    p.add_argument(
        "--iterations",
        type=int,
        default=30,
        help="Richardson-Lucy iterations (default: 30, matching the MATLAB pipeline).",
    )
    p.add_argument(
        "--engine",
        choices=("scipy", "torch", "dl2"),
        default="dl2",
        help="Deconvolution engine. 'dl2' (default) calls the "
        "DeconvolutionLab2 Java plugin via Fiji and matches the "
        "original MATLAB pipeline; requires Fiji + the DL2 plugin + a "
        "display (doesn't work headless). 'scipy' is pure-Python CPU-based "
        "Richardson-Lucy — no Fiji, no JVM required. 'torch' uses the same "
        "algorithm via PyTorch for MPS/CUDA acceleration.",
    )
    p.add_argument(
        "--fiji-dir",
        type=Path,
        default=None,
        help="Path to your Fiji.app directory (needed for --engine dl2). "
        "Overrides $FIJI_DIR.",
    )
    p.add_argument(
        "--torch-device",
        choices=("cpu", "mps", "cuda"),
        default="mps",
        help="Torch device when --engine=torch (default: mps for Apple Silicon).",
    )
    p.add_argument(
        "--no-crop-psf",
        action="store_true",
        help="Disable automatic PSF cropping. By default the PSF is "
        "cropped to its signal bounding box (+ margin) before "
        "deconvolution, which dramatically speeds up large PSFs.",
    )
    p.add_argument(
        "--psf-crop-margin",
        type=int,
        default=30,
        help="Pixel margin around the PSF signal when auto-cropping "
        "(default: 30). Ignored if --no-crop-psf is set.",
    )
    p.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=1,
        help="Number of images to deconvolve concurrently (default: 1, "
        "sequential). For --engine dl2, each job launches its own Fiji "
        "process (JVM + GUI window) — 2-4 is often a reasonable ceiling "
        "depending on CPU cores and RAM. For --engine scipy/torch this "
        "uses threads, which mainly helps if the computation isn't "
        "already saturating all CPU cores or the GPU.",
    )
    p.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable debug logging.",
    )
    return p


def deconvolve_folder(
    input_dir: Path,
    psf_path: Path,
    output_dir: Path,
    config: DeconvolutionConfig,
    max_workers: int = 1,
) -> list[Path]:
    """Deconvolve every ``.tif`` stack directly inside ``input_dir``.

    Unlike ``preprocess.deconvolve.deconvolve_channel`` (which expects the
    pipeline's ``<sample>/<channel_group>/`` layout), this just walks a
    single flat folder — for ad hoc single-folder, single-channel runs.

    ``max_workers`` controls how many images are deconvolved concurrently.
    Each runs in its own worker thread; for the ``dl2`` engine that means
    each concurrently launches its own Fiji subprocess, so wall-clock time
    for a large batch drops roughly in proportion to ``max_workers``
    (bounded by CPU/RAM). ``max_workers=1`` (the default) runs
    sequentially, matching the original behavior.
    """
    if not input_dir.is_dir():
        raise NotADirectoryError(input_dir)
    if not psf_path.is_file():
        raise FileNotFoundError(psf_path)

    output_dir.mkdir(parents=True, exist_ok=True)

    jobs: list[tuple[Path, Path]] = []
    for stack_path in sorted(input_dir.glob("*.tif")):
        if stack_path.name.endswith("_decon.tif"):
            continue
        jobs.append((stack_path, output_dir / f"{stack_path.stem}_decon.tif"))

    if not jobs:
        log.info("No .tif stacks found in %s", input_dir)
        return []

    log.info("Found %d stack(s) to deconvolve in %s (engine=%s)",
              len(jobs), input_dir, config.engine)

    # --- PSF auto-crop (benefits all engines) ---
    effective_psf_path = psf_path
    cropped_psf_tmp: Path | None = None
    cropped_psf_arr = None

    if config.crop_psf:
        raw_psf = _load_tiff_stack(psf_path)
        cropped = _autocrop_psf(raw_psf, margin=config.psf_crop_margin)
        if cropped.shape != raw_psf.shape:
            if config.engine == "dl2":
                cropped_psf_tmp = _save_cropped_psf(cropped, suffix="_folder")
                effective_psf_path = cropped_psf_tmp
            else:
                cropped_psf_arr = cropped
        else:
            if config.engine != "dl2":
                cropped_psf_arr = raw_psf

    jobs_n = max(1, int(max_workers))

    def _run_one(input_path: Path, output_path: Path) -> Path:
        if config.engine == "dl2":
            from preprocess.deconvolve_dl2 import deconvolve_file as _dl2_deconvolve_file
            log.info("Deconvolving (DL2) %s", input_path)
            _dl2_deconvolve_file(
                input_path, output_path, effective_psf_path,
                num_iter=config.num_iter,
                fiji_dir=config.fiji_dir,
            )
        else:
            log.info("Deconvolving (%s) %s", config.engine, input_path)
            _deconvolve_file_python(input_path, output_path, python_psf, config)
        log.info("Wrote %s", output_path)
        return output_path

    outputs: list[Path] = []
    try:
        python_psf = None
        if config.engine != "dl2":
            python_psf = cropped_psf_arr if cropped_psf_arr is not None else _load_tiff_stack(psf_path)
            log.info("Loaded PSF %s with shape %s", psf_path, python_psf.shape)

        if jobs_n == 1:
            for input_path, output_path in jobs:
                outputs.append(_run_one(input_path, output_path))
        else:
            log.info("Running %d job(s) with up to %d concurrent worker(s)",
                      len(jobs), jobs_n)
            with ThreadPoolExecutor(max_workers=jobs_n) as pool:
                futures = {
                    pool.submit(_run_one, input_path, output_path): output_path
                    for input_path, output_path in jobs
                }
                for future in as_completed(futures):
                    future.result()
            outputs = [output_path for _, output_path in jobs]
    finally:
        if cropped_psf_tmp is not None:
            try:
                cropped_psf_tmp.unlink()
            except OSError:
                pass

    return outputs


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    input_dir: Path = args.input_dir.resolve()
    if not input_dir.is_dir():
        parser.error(f"--input-dir is not a directory: {input_dir}")

    psf_path: Path = args.psf.resolve()
    if not psf_path.is_file():
        parser.error(f"--psf not found: {psf_path}")

    output_dir: Path = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else input_dir.parent / f"{input_dir.name}_decon"
    )

    config = DeconvolutionConfig(
        num_iter=args.iterations,
        engine=args.engine,
        torch_device=args.torch_device,
        fiji_dir=args.fiji_dir,
        crop_psf=not args.no_crop_psf,
        psf_crop_margin=args.psf_crop_margin,
    )

    if args.jobs < 1:
        parser.error("--jobs must be >= 1")

    engine_label = config.engine
    if config.engine == "torch":
        engine_label = f"torch/{config.torch_device}"

    print("\n=== Deconvolve folder ===")
    print(f"Input folder:  {input_dir}")
    print(f"PSF:           {psf_path}")
    print(f"Output folder: {output_dir}")
    print(f"Engine:        {engine_label} ({config.num_iter} iterations)")
    print(f"Jobs:          {args.jobs}\n")

    outputs = deconvolve_folder(input_dir, psf_path, output_dir, config, max_workers=args.jobs)

    print(f"\n=== Done: wrote {len(outputs)} deconvoluted stack(s) to {output_dir} ===")


if __name__ == "__main__":
    main()
