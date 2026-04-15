"""CLI entry point: prepare widefield images and run deconvolution.

Runs all preprocessing steps end-to-end:

    python preprocess.py \\
        --input-dir /path/to/main \\
        --gfp-psf /path/to/gfp_psf.tif \\
        --cy-psf  /path/to/cy_psf.tif

Individual stages can be skipped with the ``--skip-*`` flags (useful for
re-running just the deconvolution after tweaking PSFs, for example).

Input layout expected::

    <input-dir>/
        sampleA/
            gfp1/           ← single-slice TIFFs from the scope
            cy1/
            rfp1/           (optional; ignored unless stacked separately)
        sampleB/
            gfp2/
            cy2/
        ...

Output layout produced (matches the original MATLAB pipeline)::

    <input-dir>/
        Decon/
            sampleA/
                gfp/gfp1.tif      ← stacked, GFP XY-shifted
                cy/cy1.tif        ← stacked, no shift
            sampleB/...
            Deconvoluted/
                sampleA/
                    gfp/gfp1_decon.tif
                    cy/cy1_decon.tif
                sampleB/...

Point ``roi_drawing/`` at each ``Deconvoluted/<sample>/`` folder to run
ROI detection on that sample's paired gfp/cy stacks.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Make the src/ layout importable without requiring `pip install -e .`.
_SRC = Path(__file__).parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from preprocess import (  # noqa: E402
    consolidate_decon_folders,
    delete_scan_protocol_files,
    deconvolve_channel,
    rename_channel_folders,
    stack_channel_folders,
)
from preprocess.deconvolve import DeconvolutionConfig  # noqa: E402
from preprocess.stacking import DEFAULT_GFP_OFFSET  # noqa: E402


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="preprocess",
        description="Stack microscope TIFF slices and run Richardson-Lucy "
        "deconvolution (Python port of the MATLAB preprocessing pipeline).",
    )
    p.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Top-level folder containing per-sample subdirectories.",
    )
    p.add_argument(
        "--gfp-psf",
        type=Path,
        help="Measured PSF TIFF for the GFP channel. Required unless "
        "--skip-deconvolution is set.",
    )
    p.add_argument(
        "--cy-psf",
        type=Path,
        help="Measured PSF TIFF for the Cy channel. Required unless "
        "--skip-deconvolution is set.",
    )
    p.add_argument(
        "--iterations",
        type=int,
        default=30,
        help="Richardson-Lucy iterations (default: 30, matching the MATLAB pipeline).",
    )
    p.add_argument(
        "--gfp-offset",
        type=int,
        nargs=2,
        metavar=("X", "Y"),
        default=list(DEFAULT_GFP_OFFSET),
        help=f"XY pixel offset applied to every GFP slice "
        f"(default: {DEFAULT_GFP_OFFSET[0]} {DEFAULT_GFP_OFFSET[1]}). "
        "Use '0 0' to disable.",
    )
    p.add_argument(
        "--engine",
        choices=("scipy", "torch", "dl2"),
        default="scipy",
        help="Deconvolution engine. 'scipy' (default) is pure-Python CPU-based "
        "Richardson-Lucy — no Fiji, no JVM required. 'torch' uses the same "
        "algorithm via PyTorch for MPS/CUDA acceleration. 'dl2' calls the "
        "DeconvolutionLab2 Java plugin via PyImageJ and requires Fiji with "
        "the DL2 plugin installed; note that DL2 does not work in fully "
        "headless environments due to AWT GUI requirements.",
    )
    p.add_argument(
        "--fiji-dir",
        type=Path,
        default=None,
        help="Path to your Fiji.app directory (needed for --engine dl2). "
        "Overrides $FIJI_DIR. See README for install instructions.",
    )
    p.add_argument(
        "--torch-device",
        choices=("cpu", "mps", "cuda"),
        default="mps",
        help="Torch device when --engine=torch (default: mps for Apple Silicon).",
    )
    p.add_argument(
        "--sequential-channels",
        action="store_true",
        help="Deconvolve GFP then Cy sequentially. By default they run "
        "concurrently (one thread or process per channel).",
    )
    p.add_argument(
        "--skip-cleanup",
        action="store_true",
        help="Skip folder renaming and .scanprotocol deletion.",
    )
    p.add_argument(
        "--skip-stacking",
        action="store_true",
        help="Skip stacking and per-sample Decon folder creation.",
    )
    p.add_argument(
        "--skip-consolidate",
        action="store_true",
        help="Skip moving per-sample Decon folders into <input-dir>/Decon/.",
    )
    p.add_argument(
        "--skip-deconvolution",
        action="store_true",
        help="Skip Richardson-Lucy deconvolution.",
    )
    p.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable debug logging.",
    )
    return p


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    main_folder: Path = args.input_dir.resolve()
    if not main_folder.is_dir():
        parser.error(f"--input-dir is not a directory: {main_folder}")

    if not args.skip_deconvolution:
        if args.gfp_psf is None or args.cy_psf is None:
            parser.error(
                "--gfp-psf and --cy-psf are required unless --skip-deconvolution is set."
            )
        for psf_path, label in ((args.gfp_psf, "GFP"), (args.cy_psf, "Cy")):
            if not psf_path.is_file():
                parser.error(f"{label} PSF not found: {psf_path}")

    xy_offset = (int(args.gfp_offset[0]), int(args.gfp_offset[1]))

    print(f"\n=== Preprocessing pipeline ===")
    print(f"Input folder: {main_folder}")
    print()

    # 1. Cleanup
    if not args.skip_cleanup:
        print(">>> Step 1: Cleaning up folders & deleting .scanprotocol files")
        renamed = rename_channel_folders(main_folder)
        deleted = delete_scan_protocol_files(main_folder)
        print(f"    Renamed {renamed} channel folder(s). Deleted {deleted} "
              f".scanprotocol file(s).\n")
    else:
        print(">>> Step 1: SKIPPED (--skip-cleanup)\n")

    # 2. Stacking + XY offset
    if not args.skip_stacking:
        print(">>> Step 2: Stacking GFP slices with XY offset "
              f"{xy_offset} -> <sample>/Decon/gfp/")
        gfp_stacks = stack_channel_folders(
            main_folder, channel_prefix="gfp", decon_group="gfp",
            xy_offset=xy_offset if xy_offset != (0, 0) else None,
        )
        print(f"    Wrote {len(gfp_stacks)} GFP stack(s).")

        print(">>> Step 3: Stacking Cy slices (no offset) -> <sample>/Decon/cy/")
        cy_stacks = stack_channel_folders(
            main_folder, channel_prefix="cy", decon_group="cy", xy_offset=None,
        )
        print(f"    Wrote {len(cy_stacks)} Cy stack(s).\n")
    else:
        print(">>> Steps 2 & 3: SKIPPED (--skip-stacking)\n")

    # 3. Consolidate
    if not args.skip_consolidate:
        print(f">>> Step 4: Consolidating per-sample Decon/ into "
              f"{main_folder / 'Decon'}/")
        moved = consolidate_decon_folders(main_folder)
        print(f"    Consolidated {len(moved)} sample folder(s).\n")
    else:
        print(">>> Step 4: SKIPPED (--skip-consolidate)\n")

    # 4. Deconvolution
    if not args.skip_deconvolution:
        config = DeconvolutionConfig(
            num_iter=args.iterations,
            engine=args.engine,
            torch_device=args.torch_device,
            fiji_dir=args.fiji_dir,
        )
        decon_folder = main_folder / "Decon"
        if not decon_folder.is_dir():
            parser.error(
                f"Expected {decon_folder} to exist for deconvolution. "
                "Did you skip stacking or consolidate?"
            )

        engine_label = config.engine
        if config.engine == "torch":
            engine_label = f"torch/{config.torch_device}"
        print(f">>> Step 5: Deconvolution (Richardson-Lucy, "
              f"{config.num_iter} iter, engine={engine_label})")
        print(f"    GFP PSF: {args.gfp_psf}")
        print(f"    Cy  PSF: {args.cy_psf}")

        # For the DL2 engine: pre-initialize the JVM on the main thread so
        # the two channel workers share a single Fiji gateway. (JVM startup
        # isn't thread-safe and can't be restarted; one gateway per
        # process, used from multiple Java-backed threads.)
        if config.engine == "dl2":
            from preprocess.deconvolve_dl2 import _get_gateway  # noqa: WPS450
            _get_gateway(config.fiji_dir)

        # Parallelize GFP + Cy5 unless the user asked for sequential, or
        # the engine is torch on a single GPU (where they'd just serialize
        # on the device anyway).
        parallel = not args.sequential_channels and config.engine != "torch"

        if parallel:
            from concurrent.futures import ThreadPoolExecutor
            print(f"    Running GFP + Cy concurrently (engine={config.engine}).")
            with ThreadPoolExecutor(max_workers=2) as pool:
                gfp_fut = pool.submit(
                    deconvolve_channel, decon_folder, args.gfp_psf, "gfp", config,
                )
                cy_fut = pool.submit(
                    deconvolve_channel, decon_folder, args.cy_psf, "cy", config,
                )
                gfp_out = gfp_fut.result()
                cy_out = cy_fut.result()
        else:
            gfp_out = deconvolve_channel(
                decon_folder, args.gfp_psf, "gfp", config,
            )
            cy_out = deconvolve_channel(
                decon_folder, args.cy_psf, "cy", config,
            )

        print(f"    Wrote {len(gfp_out)} GFP deconvoluted stack(s).")
        print(f"    Wrote {len(cy_out)} Cy deconvoluted stack(s).\n")
    else:
        print(">>> Step 5: SKIPPED (--skip-deconvolution)\n")

    print("=== Done ===")
    if not args.skip_deconvolution:
        print(f"Deconvoluted outputs are under: "
              f"{main_folder / 'Decon' / 'Deconvoluted'}/")
        print("Point roi_drawing/ at that folder next.")


if __name__ == "__main__":
    main()
