"""Entry point: interactive-by-default CLI for cell ROI detection.

Usage (interactive, the common case):

    python roi_detect.py

Non-interactive batch mode:

    python roi_detect.py --batch path/to/folder --pan-channel gfp

Non-interactive single-pair:

    python roi_detect.py --gfp a.tif --cy5 b.tif --pan-channel gfp

Any parameter you don't specify on the command line falls back to an
interactive prompt, so you can mix and match.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import asdict
from pathlib import Path

# Make the `src/` layout importable without requiring `pip install -e .`.
_SRC = Path(__file__).parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from roi_detect import filter as roi_filter  # noqa: E402
from roi_detect import interactive, io  # noqa: E402
from roi_detect.segment import CellposeSegmenter, SegmentationParams  # noqa: E402


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="roi_detect",
        description="Detect cells that are positive in both GFP and Cy5 and "
        "export ImageJ ROIs around their full membrane.",
    )
    p.add_argument(
        "--batch",
        type=Path,
        help="Folder containing GFP/Cy5 image pairs to process in batch.",
    )
    p.add_argument("--gfp", type=Path, help="Path to a GFP TIFF (single-pair mode).")
    p.add_argument("--cy5", type=Path, help="Path to a Cy5 TIFF (single-pair mode).")
    p.add_argument(
        "--pan-channel",
        choices=("gfp", "cy5"),
        help="Which channel marks ALL cells.",
    )
    p.add_argument(
        "--diameter",
        type=float,
        help="Approximate cell diameter in pixels (omit to auto-estimate).",
    )
    p.add_argument("--blur-sigma", type=float, help="Gaussian blur sigma (default 2.0).")
    p.add_argument(
        "--gfp-threshold",
        help="Threshold spec for GFP: 'auto' (Otsu), 'pNN' (percentile), "
        "'zNN' (robust z-score, e.g. 'z2'), 'all' (skip this channel's "
        "filter), or a number.",
    )
    p.add_argument(
        "--cy5-threshold",
        help="Threshold spec for Cy5: 'auto' (Otsu), 'pNN' (percentile), "
        "'zNN' (robust z-score, e.g. 'z2'), 'all' (skip this channel's "
        "filter), or a number.",
    )
    p.add_argument(
        "--flow-threshold",
        type=float,
        help="Cellpose flow_threshold (default 0.4; higher = more irregular shapes allowed).",
    )
    p.add_argument(
        "--cellprob-threshold",
        type=float,
        help="Cellpose cellprob_threshold (default 0.0; lower = catches dimmer cell edges).",
    )
    p.add_argument(
        "--niter",
        type=int,
        help="Cellpose dynamics iterations (default auto; use 2000+ for long thin processes).",
    )
    p.add_argument(
        "--min-size",
        type=int,
        help="Minimum detected object size in pixels (default 15).",
    )
    p.add_argument(
        "--segment-on",
        choices=("pan", "both"),
        help="'pan' = segment on pan-cell channel only; 'both' = segment on max "
        "of both channels so ROIs cover signal in either channel.",
    )
    p.add_argument(
        "--intensity-metric",
        help="Which per-cell intensity to threshold on: 'mean' (default), "
        "'topN' (e.g. 'top95'), 'brightN' (e.g. 'bright10'), or "
        "'contrast' (ratio of bright pixels to median within each cell — "
        "best for punctate signal with variable autofluorescence).",
    )
    p.add_argument(
        "--irregular-shapes",
        action="store_true",
        help="Preset for non-round cells (INS-1, neurons, cells with processes).",
    )
    p.add_argument(
        "--output-dir", type=Path, help="Where to write results (default: input folder)."
    )
    p.add_argument(
        "--no-gpu",
        action="store_true",
        help="Force CPU mode even if a GPU is available.",
    )
    p.add_argument(
        "--non-interactive",
        action="store_true",
        help="Fail instead of prompting if any required parameter is missing.",
    )
    return p


def merge_cli_into_settings(
    args: argparse.Namespace, settings: interactive.RunSettings
) -> interactive.RunSettings:
    # Apply the irregular-shapes preset first so explicit CLI flags below
    # can still override individual values from the preset.
    if args.irregular_shapes:
        from roi_detect.segment import irregular_shapes_preset, SegmentationParams

        preset = irregular_shapes_preset(SegmentationParams())
        settings.irregular_shapes = True
        settings.segment_on = "both"
        settings.blur_sigma = preset.blur_sigma
        settings.flow_threshold = preset.flow_threshold
        settings.cellprob_threshold = preset.cellprob_threshold
        settings.niter = preset.niter

    if args.pan_channel:
        settings.pan_channel = args.pan_channel
    if args.diameter is not None:
        settings.diameter = args.diameter
    if args.blur_sigma is not None:
        settings.blur_sigma = args.blur_sigma
    if args.gfp_threshold is not None:
        settings.gfp_threshold = args.gfp_threshold
    if args.cy5_threshold is not None:
        settings.cy5_threshold = args.cy5_threshold
    if args.flow_threshold is not None:
        settings.flow_threshold = args.flow_threshold
    if args.cellprob_threshold is not None:
        settings.cellprob_threshold = args.cellprob_threshold
    if args.niter is not None:
        settings.niter = args.niter
    if args.min_size is not None:
        settings.min_size = args.min_size
    if args.segment_on is not None:
        settings.segment_on = args.segment_on
    if args.intensity_metric is not None:
        settings.intensity_metric = args.intensity_metric
    if args.no_gpu:
        settings.use_gpu = False
    return settings


def process_pair(
    pair: io.ImagePair,
    segmenter: CellposeSegmenter,
    settings: interactive.RunSettings,
    output_dir: Path,
    pair_index: int | None = None,
    roi_dir: Path | None = None,
) -> dict:
    """Run the full pipeline on one image pair. Returns a summary dict.

    If ``pair_index`` is given, the ROI zip is also saved as
    ``<roi_dir>/roi_original/<NN>.zip``. ``roi_dir`` defaults to
    ``output_dir`` if not specified.
    """
    from roi_detect.segment import preprocess_composite

    gfp, cy5 = io.load_pair(pair)

    if settings.segment_on == "both":
        composite = preprocess_composite(gfp, cy5, blur_sigma=segmenter.params.blur_sigma)
        masks = segmenter.segment(composite, already_preprocessed=True)
    else:
        pan = gfp if settings.pan_channel == "gfp" else cy5
        masks = segmenter.segment(pan)

    # Parse intensity metric to get the parameters for measure_cells.
    metric_kind, metric_pct = roi_filter._parse_metric(settings.intensity_metric)
    top_pct = metric_pct if metric_kind == "top" else 95.0
    bright_pct = metric_pct if metric_kind in ("bright", "contrast") else 10.0

    rows, _gfp_bg, _cy5_bg = roi_filter.measure_cells(
        masks, gfp, cy5, top_percentile=top_pct, bright_pct=bright_pct,
    )
    measurements, gfp_thr, cy5_thr = roi_filter.score_cells(
        rows, settings.gfp_threshold, settings.cy5_threshold,
        intensity_metric=settings.intensity_metric,
    )

    kept_labels = {m.label for m in measurements if m.kept}

    # Build a kept-only mask for the saved mask TIFF.
    import numpy as np

    kept_mask = np.zeros_like(masks)
    for label in kept_labels:
        kept_mask[masks == label] = label

    stem = pair.stem
    output_dir.mkdir(parents=True, exist_ok=True)
    roi_path = output_dir / f"{stem}_rois.zip"
    mask_path = output_dir / f"{stem}_masks.tif"
    overlay_path = output_dir / f"{stem}_overlay.png"
    csv_path = output_dir / f"{stem}_measurements.csv"

    io.save_roi_zip(roi_path, masks, sorted(kept_labels))
    io.save_masks(mask_path, kept_mask)
    io.save_overlay(overlay_path, gfp, cy5, masks, kept_labels)

    # Save a numbered copy in the shared roi_original/ folder. This is the
    # untouched automated output; the user edits it in the ImageJ
    # Color_Merge macro and writes the edited version to a sibling roi/
    # folder, which the next pipeline stage (roi_cropping) reads from.
    if pair_index is not None:
        roi_base = roi_dir if roi_dir is not None else output_dir
        roi_folder = roi_base / "roi_original"
        roi_folder.mkdir(parents=True, exist_ok=True)
        numbered_path = roi_folder / f"{pair_index:02d}.zip"
        io.save_roi_zip(numbered_path, masks, sorted(kept_labels))

    csv_rows = []
    for m in measurements:
        row = asdict(m)
        row["gfp_threshold"] = gfp_thr
        row["cy5_threshold"] = cy5_thr
        csv_rows.append(row)
    io.save_measurements_csv(csv_path, csv_rows)

    return {
        "pair": pair.describe(),
        "total_cells": len(measurements),
        "kept_cells": len(kept_labels),
        "gfp_threshold": gfp_thr,
        "cy5_threshold": cy5_thr,
    }


def run_batch(
    folder: Path,
    settings: interactive.RunSettings,
    output_dir: Path | None,
    interactive_mode: bool,
) -> None:
    pairs, detected = io.find_pairs(folder)
    if not pairs:
        print(
            f"No image pairs found in {folder}. "
            "Expected TIFFs named e.g. '<stem>_GFP.tif' + '<stem>_Cy5.tif', "
            "or 2-channel multi-page TIFFs."
        )
        return

    if detected:
        print(f"\nDetected naming convention: {detected}")
    print(f"Found {len(pairs)} image pair(s) in {folder}:")
    for pair in pairs:
        print(f"  {pair.describe()}")

    # Handle multi-channel files: ask once which index is GFP, apply to all.
    multi_pairs = [p for p in pairs if p.multi_path is not None]
    if multi_pairs:
        print(
            f"\n{len(multi_pairs)} file(s) are multi-channel and need channel "
            "assignment."
        )
        if interactive_mode:
            gfp_idx = interactive.ask_int("GFP channel index", 0) or 0
            cy5_idx = interactive.ask_int("Cy5 channel index", 1) or 1
        else:
            gfp_idx, cy5_idx = 0, 1
        for p in multi_pairs:
            p.gfp_index = gfp_idx
            p.cy5_index = cy5_idx

    if interactive_mode and not interactive.ask_yes_no("\nProceed?", default=True):
        print("Aborted.")
        return

    if output_dir is None:
        output_dir = folder / "roi_results"
    print(f"\nWriting results to: {output_dir}")
    print(f"ROI zips will be in: {folder / 'roi_original'}")

    try:
        from tqdm import tqdm
    except ImportError:
        def tqdm(x, **_kwargs):
            return x

    segmenter = CellposeSegmenter(_make_segmentation_params(settings))

    summaries = []
    for idx, pair in enumerate(tqdm(pairs, desc="Segmenting", unit="pair"), start=1):
        try:
            summary = process_pair(
                pair, segmenter, settings, output_dir,
                pair_index=idx, roi_dir=folder,
            )
            summaries.append(summary)
        except Exception as e:
            print(f"\n! Failed on {pair.describe()}: {e}")

    print("\n--- Batch summary ---")
    total_cells = sum(s["total_cells"] for s in summaries)
    total_kept = sum(s["kept_cells"] for s in summaries)
    print(
        f"{len(summaries)}/{len(pairs)} pair(s) processed. "
        f"{total_kept}/{total_cells} cells kept as double-positive."
    )
    for s in summaries:
        print(f"  {s['pair']}: {s['kept_cells']}/{s['total_cells']} kept")


def run_single(
    pair: io.ImagePair,
    settings: interactive.RunSettings,
    output_dir: Path | None,
) -> None:
    anchor = pair.multi_path or pair.gfp_path
    input_folder = anchor.parent  # type: ignore[union-attr]
    if output_dir is None:
        output_dir = input_folder / "roi_results"

    print(f"\nProcessing: {pair.describe()}")
    print(f"Writing results to: {output_dir}")

    segmenter = CellposeSegmenter(_make_segmentation_params(settings))
    summary = process_pair(
        pair, segmenter, settings, output_dir,
        pair_index=1, roi_dir=input_folder,
    )
    print(
        f"\nDone: {summary['kept_cells']}/{summary['total_cells']} cells kept "
        f"(GFP threshold={summary['gfp_threshold']:.3f}, "
        f"Cy5 threshold={summary['cy5_threshold']:.3f})"
    )


def _make_segmentation_params(settings: interactive.RunSettings) -> SegmentationParams:
    return SegmentationParams(
        diameter=settings.diameter,
        blur_sigma=settings.blur_sigma,
        flow_threshold=settings.flow_threshold,
        cellprob_threshold=settings.cellprob_threshold,
        niter=settings.niter,
        min_size=settings.min_size,
        use_gpu=settings.use_gpu,
    )


def _find_partner_file(file: Path) -> io.ImagePair | None:
    """Try to find the matching partner file in the same folder as ``file``."""
    pairs, _ = io.find_pairs(file.parent)
    for pair in pairs:
        if pair.multi_path == file:
            return pair
        if pair.gfp_path == file or pair.cy5_path == file:
            return pair
    return None


def interactive_mode_main(
    args: argparse.Namespace, settings: interactive.RunSettings
) -> None:
    print("=" * 60)
    print(" ROI Detection — cells positive in GFP + Cy5")
    print("=" * 60)

    print("\nMode?")
    print("  [1] Batch: process a whole folder of image pairs")
    print("  [2] Single: process one image pair (for testing)")
    choice = interactive.ask("Choice", default="1")

    if choice == "2":
        initial = settings.last_folder or ""
        file = interactive.pick_file("Pick any image in the pair", initial)
        if file is None:
            print("No file selected. Exiting.")
            return

        pair = _find_partner_file(file)
        if pair is None:
            print(
                f"Couldn't auto-find a partner for {file.name}. "
                "Please pick the second file."
            )
            second = interactive.pick_file("Pick the second image", str(file.parent))
            if second is None:
                print("No file selected. Exiting.")
                return
            pair = io.ImagePair(gfp_path=file, cy5_path=second)

        if pair.multi_path is not None and (
            pair.gfp_index is None or pair.cy5_index is None
        ):
            pair.gfp_index = interactive.ask_int("GFP channel index", 0) or 0
            pair.cy5_index = interactive.ask_int("Cy5 channel index", 1) or 1

        settings.last_folder = str(file.parent)
        settings = interactive.prompt_run_settings(settings)
        settings.save()
        run_single(pair, settings, args.output_dir)
        return

    # Default: batch mode.
    initial = settings.last_folder or ""
    folder = interactive.pick_folder("Pick a folder of image pairs", initial)
    if folder is None:
        print("No folder selected. Exiting.")
        return
    settings.last_folder = str(folder)

    settings = interactive.prompt_run_settings(settings)
    settings.save()
    run_batch(folder, settings, args.output_dir, interactive_mode=True)


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    settings = interactive.RunSettings.load()
    settings = merge_cli_into_settings(args, settings)

    # Decide which branch to take.
    if args.batch is not None:
        if not args.pan_channel and not args.non_interactive:
            settings = interactive.prompt_run_settings(settings)
        elif not args.pan_channel:
            parser.error("--pan-channel is required in non-interactive batch mode")
        settings.save()
        run_batch(args.batch, settings, args.output_dir, interactive_mode=not args.non_interactive)
        return

    if args.gfp is not None and args.cy5 is not None:
        if not args.pan_channel and not args.non_interactive:
            settings = interactive.prompt_run_settings(settings)
        elif not args.pan_channel:
            parser.error("--pan-channel is required in non-interactive single-pair mode")
        settings.save()
        pair = io.ImagePair(gfp_path=args.gfp, cy5_path=args.cy5)
        run_single(pair, settings, args.output_dir)
        return

    if args.non_interactive:
        parser.error(
            "No input specified. Pass --batch or --gfp/--cy5, "
            "or drop --non-interactive to be prompted."
        )

    interactive_mode_main(args, settings)


if __name__ == "__main__":
    main()
