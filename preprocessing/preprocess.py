"""CLI entry point: prepare widefield images and run deconvolution.

Runs all preprocessing steps end-to-end:

    python preprocess.py \\
        --input-dir /path/to/main \\
        --gfp-psf /path/to/gfp_psf.tif \\
        --cy-psf  /path/to/cy_psf.tif

Individual stages can be skipped with the ``--skip-*`` flags (useful for
re-running just the deconvolution after tweaking PSFs, for example).

Pass ``--jobs N`` to deconvolve multiple images within each channel
concurrently (default: 1, sequential) — the same concurrency model as
``deconvolve_folder.py --jobs``. Combines multiplicatively with
channel-level concurrency (channels run concurrently by default too; pass
``--sequential-channels`` to disable that separately).

Input layout expected::

    <input-dir>/
        sampleA/
            gfp1/           ← single-slice TIFFs from the scope
            cy1/            (GFP/Cy experiments)
            rfp1/           (GFP/RFP experiments)
        sampleB/
            gfp2/
            cy2/
        ...

A sample normally pairs GFP with Cy *or* RFP, not both — pass ``--cy-psf``
and/or ``--rfp-psf`` depending on which second channel this dataset uses.
Both are stacked opportunistically; only channels with a PSF supplied are
deconvolved.

**GFP XY offset**: the chromatic registration shift depends on which
channel GFP is paired with. Default ``--gfp-offset`` (``5 -2``) is
calibrated for GFP/Cy5. For GFP/RFP, pass ``--gfp-offset -3 1`` explicitly.

Output layout produced (matches the original MATLAB pipeline)::

    <input-dir>/
        Decon/
            sampleA/
                gfp/gfp1.tif      ← stacked, GFP XY-shifted
                cy/cy1.tif        ← stacked, no shift
                rfp/rfp1.tif      ← stacked, no shift (if present)
            sampleB/...
            Deconvoluted/
                sampleA/
                    gfp/gfp1_decon.tif
                    cy/cy1_decon.tif
                    rfp/rfp1_decon.tif
                sampleB/...

Point ``roi_drawing/`` at each ``Deconvoluted/<sample>/`` folder to run
ROI detection on that sample's paired gfp/cy (or gfp/rfp) stacks.
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
    zero_pad_channel_folders,
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
        help="Measured PSF TIFF for the Cy channel. At least one of "
        "--cy-psf / --rfp-psf is required unless --skip-deconvolution is set.",
    )
    p.add_argument(
        "--rfp-psf",
        type=Path,
        help="Measured PSF TIFF for the RFP channel. At least one of "
        "--cy-psf / --rfp-psf is required unless --skip-deconvolution is set.",
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
        f"(default: {DEFAULT_GFP_OFFSET[0]} {DEFAULT_GFP_OFFSET[1]}, "
        "calibrated for GFP/Cy5). For GFP/RFP pass '-3 1' instead. "
        "Use '0 0' to disable.",
    )
    p.add_argument(
        "--engine",
        choices=("scipy", "torch", "dl2"),
        default="dl2",
        help="Deconvolution engine. 'dl2' (default) calls the "
        "DeconvolutionLab2 Java plugin via PyImageJ and matches the "
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
        "Overrides $FIJI_DIR. See README for install instructions.",
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
        "--sequential-channels",
        action="store_true",
        help="Deconvolve channels one at a time instead of concurrently "
        "(one thread per channel by default).",
    )
    p.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=1,
        help="Number of images WITHIN each channel to deconvolve "
        "concurrently (default: 1, sequential) — same concurrency model as "
        "deconvolve_folder.py --jobs. For --engine dl2 each job launches "
        "its own Fiji process; 2-4 is often a reasonable ceiling depending "
        "on CPU cores and RAM. This is independent of channel-level "
        "concurrency (--sequential-channels): with the default settings "
        "and --jobs 3 on a GFP/RFP run, up to 2 channels x 3 jobs = 6 Fiji "
        "processes can run at once, so lower --jobs or pass "
        "--sequential-channels on memory-constrained machines.",
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
        if args.gfp_psf is None:
            parser.error("--gfp-psf is required unless --skip-deconvolution is set.")
        if args.cy_psf is None and args.rfp_psf is None:
            parser.error(
                "At least one of --cy-psf / --rfp-psf is required unless "
                "--skip-deconvolution is set."
            )
        psfs_to_check = [(args.gfp_psf, "GFP")]
        if args.cy_psf is not None:
            psfs_to_check.append((args.cy_psf, "Cy"))
        if args.rfp_psf is not None:
            psfs_to_check.append((args.rfp_psf, "RFP"))
        for psf_path, label in psfs_to_check:
            if not psf_path.is_file():
                parser.error(f"{label} PSF not found: {psf_path}")

    if args.jobs < 1:
        parser.error("--jobs must be >= 1")

    xy_offset = (int(args.gfp_offset[0]), int(args.gfp_offset[1]))

    print(f"\n=== Preprocessing pipeline ===")
    print(f"Input folder: {main_folder}")
    print()

    # 1. Cleanup
    if not args.skip_cleanup:
        print(">>> Step 1: Cleaning up folders & deleting .scanprotocol files")
        renamed = rename_channel_folders(main_folder)
        padded = zero_pad_channel_folders(main_folder)
        deleted = delete_scan_protocol_files(main_folder)
        print(f"    Renamed {renamed} channel folder(s). "
              f"Zero-padded {padded} folder(s). "
              f"Deleted {deleted} .scanprotocol file(s).\n")
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

        print(">>> Step 4: Stacking RFP slices (no offset) -> <sample>/Decon/rfp/")
        rfp_stacks = stack_channel_folders(
            main_folder, channel_prefix="rfp", decon_group="rfp", xy_offset=None,
        )
        print(f"    Wrote {len(rfp_stacks)} RFP stack(s).\n")
    else:
        print(">>> Steps 2-4: SKIPPED (--skip-stacking)\n")

    # 3. Consolidate
    if not args.skip_consolidate:
        print(f">>> Step 5: Consolidating per-sample Decon/ into "
              f"{main_folder / 'Decon'}/")
        moved = consolidate_decon_folders(main_folder)
        print(f"    Consolidated {len(moved)} sample folder(s).\n")
    else:
        print(">>> Step 5: SKIPPED (--skip-consolidate)\n")

    # 4. Deconvolution
    if not args.skip_deconvolution:
        config = DeconvolutionConfig(
            num_iter=args.iterations,
            engine=args.engine,
            torch_device=args.torch_device,
            fiji_dir=args.fiji_dir,
            crop_psf=not args.no_crop_psf,
            psf_crop_margin=args.psf_crop_margin,
        )
        decon_folder = main_folder / "Decon"
        if not decon_folder.is_dir():
            parser.error(
                f"Expected {decon_folder} to exist for deconvolution. "
                "Did you skip stacking or consolidate?"
            )

        # Channels to deconvolve: GFP always, Cy/RFP whichever have a PSF.
        channels = [("GFP", "gfp", args.gfp_psf)]
        if args.cy_psf is not None:
            channels.append(("Cy", "cy", args.cy_psf))
        if args.rfp_psf is not None:
            channels.append(("RFP", "rfp", args.rfp_psf))

        engine_label = config.engine
        if config.engine == "torch":
            engine_label = f"torch/{config.torch_device}"
        print(f">>> Step 6: Deconvolution (Richardson-Lucy, "
              f"{config.num_iter} iter, engine={engine_label}, jobs={args.jobs})")
        for label, _group, psf_path in channels:
            print(f"    {label} PSF: {psf_path}")

        # Parallelize channels unless the user asked for sequential, or the
        # engine is torch on a single GPU (where they'd just serialize on
        # the device anyway). Independently, --jobs controls how many
        # images WITHIN each channel run concurrently, so the two combine
        # multiplicatively — see the --jobs help text.
        parallel = not args.sequential_channels and config.engine != "torch"

        outputs: dict[str, list[Path]] = {}
        if parallel:
            from concurrent.futures import ThreadPoolExecutor
            channel_names = " + ".join(label for label, _g, _p in channels)
            print(f"    Running {channel_names} concurrently (engine={config.engine}), "
                  f"{args.jobs} job(s) per channel.")
            with ThreadPoolExecutor(max_workers=len(channels)) as pool:
                futures = {
                    label: pool.submit(
                        deconvolve_channel, decon_folder, psf_path, group, config,
                        args.jobs,
                    )
                    for label, group, psf_path in channels
                }
                outputs = {label: fut.result() for label, fut in futures.items()}
        else:
            outputs = {
                label: deconvolve_channel(decon_folder, psf_path, group, config, args.jobs)
                for label, group, psf_path in channels
            }

        for label, _group, _psf_path in channels:
            print(f"    Wrote {len(outputs[label])} {label} deconvoluted stack(s).")
        print()
    else:
        print(">>> Step 6: SKIPPED (--skip-deconvolution)\n")

    print("=== Done ===")
    if not args.skip_deconvolution:
        print(f"Deconvoluted outputs are under: "
              f"{main_folder / 'Decon' / 'Deconvoluted'}/")
        print("Point roi_drawing/ at that folder next.")


if __name__ == "__main__":
    main()
