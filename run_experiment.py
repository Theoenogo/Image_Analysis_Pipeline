"""Top-level batch runner for the Python half of the pipeline.

Chains the three stages that share a consistent folder convention:

    roi_cropping  →  background_subtraction  →  manders_mcc

This is purely a convenience wrapper. The three individual CLIs
(``roi_cropping/roi_crop.py``, ``background_subtraction/bg_subtract.py``,
``manders_mcc/standalone_analysis.py``) still work on their own and are
unchanged — use them when you want to iterate on a single stage.

Use this script when you've finished the manual ROI-editing step in
Fiji and want to push every sample in the experiment through to final
colocalization CSVs in one command.

Usage
-----

    python run_experiment.py --input-dir /path/to/main_folder

The walker matches every stage's individual behavior:

- **roi_cropping**: every folder under ``--input-dir`` that contains
  ``gfp/``, ``cy/``, and ``roi/`` subfolders is treated as a sample.
  Writes ``<sample>/Cropped/{gfp,cy,roi}/`` + ``roi.zip``.
- **background_subtraction**: every ``Cropped/`` folder produced above
  gets a sibling ``Background_Subtracted/`` written next to it.
- **manders_mcc**: every ``Background_Subtracted/`` folder gets a pair
  of CSVs written inside it (``coloc_results.csv`` +
  ``thresholds_multislice.csv``).

Any stage can be skipped with ``--skip-*`` flags (useful for re-running
only the colocalization analysis after tweaking thresholds, for example).
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).parent.resolve()

# Make each tool's src/ layout importable without pip-installing them.
for _sub in ("roi_cropping", "background_subtraction"):
    _src = _REPO_ROOT / _sub / "src"
    if _src.is_dir() and str(_src) not in sys.path:
        sys.path.insert(0, str(_src))

# manders_mcc has no src/ layout — its module is at the folder root.
_manders_dir = _REPO_ROOT / "manders_mcc"
if _manders_dir.is_dir() and str(_manders_dir) not in sys.path:
    sys.path.insert(0, str(_manders_dir))

from roi_crop import run_pipeline as run_roi_crop  # noqa: E402
from bg_subtract import run_pipeline as run_bg_subtract  # noqa: E402
from standalone_analysis import run_pipeline as run_manders  # noqa: E402


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_experiment",
        description=(
            "Batch-run the roi_cropping → background_subtraction → "
            "manders_mcc stages across every sample under --input-dir."
        ),
    )
    p.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Top-level folder to walk. Each sample (folder with gfp/, cy/, "
        "roi/ subfolders) is processed through all three stages.",
    )

    # roi_cropping options
    p.add_argument("--gfp-dirname", default="gfp",
                   help="GFP channel folder name (default: gfp).")
    p.add_argument("--cy-dirname", default="cy",
                   help="Cy channel folder name (default: cy).")
    p.add_argument("--roi-dirname", default="roi",
                   help="Edited ROI zip folder name (default: roi).")
    p.add_argument("--gfp-prefix", default="gfp",
                   help="Cropped GFP filename prefix (default: gfp).")
    p.add_argument("--cy-prefix", default="cy",
                   help="Cropped Cy filename prefix (default: cy).")

    # background_subtraction options
    p.add_argument("--gfp-multiplier", type=float, default=1.25,
                   help="GFP per-cell mean multiplier (default: 1.25, matches "
                        "the ImageJ macro default).")
    p.add_argument("--cy-multiplier", type=float, default=1.25,
                   help="Cy per-cell mean multiplier (default: 1.25).")
    p.add_argument("--bg-floor", type=float, default=100.0,
                   help="Lower clamp on the subtraction value (default: 100).")
    p.add_argument("--bg-ceiling", type=float, default=5000.0,
                   help="Upper clamp on the subtraction value (default: 5000).")

    # Stage-skipping flags
    p.add_argument("--skip-roi-cropping", action="store_true",
                   help="Skip the roi_cropping stage (expects Cropped/ folders "
                        "to already exist).")
    p.add_argument("--skip-background-subtraction", action="store_true",
                   help="Skip the background-subtraction stage (expects "
                        "Background_Subtracted/ folders to already exist).")
    p.add_argument("--skip-manders-mcc", action="store_true",
                   help="Skip the colocalization-analysis stage.")

    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose logging.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    input_dir: Path = args.input_dir.resolve()
    if not input_dir.is_dir():
        print(f"Error: --input-dir is not a directory: {input_dir}")
        return 1

    print(f"\n=== Experiment batch run ===")
    print(f"Input folder: {input_dir}\n")

    # ---- Stage 1: roi_cropping ----
    if not args.skip_roi_cropping:
        print(">>> Stage 1: ROI cropping")
        crop_results = run_roi_crop(
            input_dir,
            gfp_dirname=args.gfp_dirname,
            cy_dirname=args.cy_dirname,
            roi_dirname=args.roi_dirname,
            gfp_prefix=args.gfp_prefix,
            cy_prefix=args.cy_prefix,
        )
        total_cells = sum(len(o.gfp) for o in crop_results.values())
        print(f"    Cropped {total_cells} cell(s) across "
              f"{len(crop_results)} sample(s).\n")
    else:
        print(">>> Stage 1: SKIPPED (--skip-roi-cropping)\n")

    # ---- Stage 2: background_subtraction ----
    if not args.skip_background_subtraction:
        print(">>> Stage 2: Background subtraction")
        bg_results = run_bg_subtract(
            input_dir,
            gfp_multiplier=args.gfp_multiplier,
            cy_multiplier=args.cy_multiplier,
            floor=args.bg_floor,
            ceiling=args.bg_ceiling,
        )
        total_pairs = sum(len(o.gfp) for o in bg_results.values())
        print(f"    Background-subtracted {total_pairs} pair(s) across "
              f"{len(bg_results)} Cropped/ folder(s).\n")
    else:
        print(">>> Stage 2: SKIPPED (--skip-background-subtraction)\n")

    # ---- Stage 3: manders_mcc ----
    if not args.skip_manders_mcc:
        print(">>> Stage 3: Manders colocalization analysis")
        manders_results = run_manders(input_dir, verbose=args.verbose)
        if not manders_results:
            print("    No Background_Subtracted/ folders found — nothing to "
                  "analyze.")
        else:
            total_pairs = sum(r["pairs_processed"] for r in manders_results.values())
            total_slices = sum(r["slices_processed"] for r in manders_results.values())
            print(f"\n    Analyzed {total_pairs} image pair(s) / "
                  f"{total_slices} slice(s) across "
                  f"{len(manders_results)} Background_Subtracted/ folder(s).")
            for folder, summary in manders_results.items():
                print(f"      {folder} → {Path(summary['results_csv']).name}")
    else:
        print(">>> Stage 3: SKIPPED (--skip-manders-mcc)\n")

    print("\n=== Done ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
