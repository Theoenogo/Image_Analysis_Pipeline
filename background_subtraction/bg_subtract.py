"""CLI entry point for background_subtraction.

Python port of ``Background_Subtraction_Tailored.ijm``. Subtracts a
per-cell constant (``mean_intensity_inside_ROI * multiplier``, clamped to
``[100, 5000]``) from every slice of every cropped stack.

Input layout (produced by ``roi_cropping/``)::

    <sample>/Cropped/
        gfp/gfp1.tif, gfp2.tif, ...
        cy/cy1.tif, cy2.tif, ...
        roi.zip

Output layout (consumed by ``manders_mcc/``)::

    <sample>/Background_Subtracted/
        gfp/gfp1.tif, gfp2.tif, ...
        cy/cy1.tif, cy2.tif, ...
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

from bg_subtract import run_pipeline  # noqa: E402
from bg_subtract.pipeline import subtract_sample_folder  # noqa: E402


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="bg_subtract",
        description=(
            "Per-cell background subtraction. Walks input-dir for "
            "Cropped/ folders and writes Background_Subtracted/ siblings."
        ),
    )
    p.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Folder to walk recursively for Cropped/ subfolders. In "
        "--direct mode, this is instead the single folder to process "
        "(must directly contain gfp/, cy/, and an ROI source).",
    )
    p.add_argument(
        "--direct",
        action="store_true",
        help="Treat --input-dir as the sample folder itself instead of "
        "walking for Cropped/ subfolders — for running background "
        "subtraction on one specific folder (e.g. output from "
        "extras/remove_body.py) that isn't named 'Cropped'.",
    )
    p.add_argument(
        "--roi-zip", default="roi.zip",
        help="ROI zip filename to look for inside each sample (default: "
        "roi.zip) — one ROI per image, matched by index.",
    )
    p.add_argument(
        "--roi-dirname", default="roi",
        help="Fallback ROI folder name if --roi-zip isn't found (default: "
        "roi) — one .roi/.zip file per image, matched positionally by "
        "filename number (e.g. extras/remove_body.py's roi/ input).",
    )
    p.add_argument("--gfp-multiplier", type=float, default=1.25,
                   help="Multiplier on the GFP per-cell mean (default: 1.25, "
                        "matches the ImageJ macro default).")
    p.add_argument("--cy-multiplier", type=float, default=1.25,
                   help="Multiplier on the CY per-cell mean (default: 1.25).")
    p.add_argument("--floor", type=float, default=100.0,
                   help="Lower clamp on the subtraction value (default: 100).")
    p.add_argument("--ceiling", type=float, default=5000.0,
                   help="Upper clamp on the subtraction value (default: 5000).")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if args.direct:
        sample_dir = args.input_dir.resolve()
        if not sample_dir.is_dir():
            print(f"--input-dir is not a directory: {sample_dir}", file=sys.stderr)
            return 1
        print(f"Direct mode: {sample_dir}")
        outs = subtract_sample_folder(
            sample_dir,
            roi_zip_name=args.roi_zip,
            roi_dirname=args.roi_dirname,
            gfp_multiplier=args.gfp_multiplier,
            cy_multiplier=args.cy_multiplier,
            floor=args.floor,
            ceiling=args.ceiling,
        )
        out_root = sample_dir.parent / "Background_Subtracted"
        print(f"\nDone. Background-subtracted {len(outs.gfp)} pair(s) -> {out_root}")
        return 0

    results = run_pipeline(
        args.input_dir,
        roi_zip_name=args.roi_zip,
        roi_dirname=args.roi_dirname,
        gfp_multiplier=args.gfp_multiplier,
        cy_multiplier=args.cy_multiplier,
        floor=args.floor,
        ceiling=args.ceiling,
    )

    total_pairs = sum(len(o.gfp) for o in results.values())
    print(f"\nDone. Background-subtracted {total_pairs} pair(s) "
          f"across {len(results)} sample(s).")
    for cropped_root, outs in results.items():
        out_root = cropped_root.parent / "Background_Subtracted"
        print(f"  {cropped_root.parent}: {len(outs.gfp)} pair(s) -> {out_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
