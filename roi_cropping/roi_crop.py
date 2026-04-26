"""CLI entry point for roi_cropping.

Replaces the MATLAB ``H_roi_duplicate_image_all_channels_Recursive.m`` plus
the ImageJ ``ROI_Select_Duplicate_TIFF_Loop.ijm`` two-step flow with a
single Python pass.

Per-sample input layout::

    <sample>/
        gfp/gfpN_decon.tif
        cy/cyN_decon.tif
        roi/NN.zip                 ← user-edited ROIs (one zip per pair)

Output layout::

    <sample>/Cropped/
        gfp/gfp1.tif, gfp2.tif, ...   ← one cropped stack per ROI
        cy/cy1.tif, cy2.tif, ...
        roi/1.roi, 2.roi, ...         ← per-cell ROI in bbox-local coords
        roi.zip                       ← combined zip in image-paired order

The ``Cropped/`` tree is what the next stage (``background_subtraction/``)
reads from.
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

from roi_crop import run_pipeline  # noqa: E402


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="roi_crop",
        description=(
            "Crop single-cell stacks from decon outputs using user-edited "
            "ROI zips. Python port of MATLAB H_roi_duplicate + ImageJ "
            "ROI_Select_Duplicate_TIFF_Loop, collapsed into one pass."
        ),
    )
    p.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help=(
            "Folder to walk. Any folder containing all of <gfp>/, <cy>/, "
            "and <roi>/ subfolders is treated as a sample."
        ),
    )
    p.add_argument("--gfp-dirname", default="gfp",
                   help="GFP channel folder name (default: gfp).")
    p.add_argument("--cy-dirname", default="cy",
                   help="Cy channel folder name (default: cy).")
    p.add_argument("--roi-dirname", default="roi",
                   help="ROI zip folder name (default: roi).")
    p.add_argument("--gfp-prefix", default="gfp",
                   help="Filename prefix for cropped GFP outputs (default: gfp).")
    p.add_argument("--cy-prefix", default="cy",
                   help="Filename prefix for cropped Cy outputs (default: cy).")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose logging.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    results = run_pipeline(
        args.input_dir,
        gfp_dirname=args.gfp_dirname,
        cy_dirname=args.cy_dirname,
        roi_dirname=args.roi_dirname,
        gfp_prefix=args.gfp_prefix,
        cy_prefix=args.cy_prefix,
    )

    total_cells = sum(len(o.gfp) for o in results.values())
    print(f"\nDone. Cropped {total_cells} cell(s) across {len(results)} sample(s).")
    for sample, outs in results.items():
        print(f"  {sample}: {len(outs.gfp)} cell(s) -> {sample / 'Cropped'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
