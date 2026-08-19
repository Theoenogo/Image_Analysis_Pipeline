"""Flat background subtraction using a low ROI percentile instead of the mean.

``bg_subtract.py`` computes its per-slice subtraction constant as
``mean_intensity_inside_ROI * multiplier``. That statistic is exactly the
one a small, very bright, spatially-concentrated region (e.g. the Golgi in
a proinsulin/insulin channel) distorts most — a handful of very bright
pixels pull the whole-ROI mean up, inflating the constant, which then gets
subtracted everywhere in the image, disproportionately erasing real, dimmer
signal far from that bright region (puncta out in cell processes / near
the plasma membrane).

This script keeps the exact same algorithm shape — one flat constant per
slice, applied to the whole image, clamped to ``[floor, ceiling]`` — and
changes exactly one thing: the constant is computed from a **low
percentile** of the ROI's pixel values instead of the mean. A low
percentile (default: 10th) reflects the typical dim background pixel
regardless of how bright a small region elsewhere in the ROI gets, since
that region only ever occupies a minority of the ROI's pixels. Same
statistic-in-a-formula shape as the method it's replacing, just a
different statistic — deliberately the smallest change that fixes the
Golgi-inflation problem, rather than a new algorithm.

(An earlier spatially-local rolling-ball approach was tried and dropped:
it correctly decouples a bright Golgi from the rest of the image, but a
rolling ball can't remove a genuinely flat, uniform background — there's
no local trend for it to roll under — so it left far-field noise almost
entirely unsuppressed. The percentile fix above solves the actual problem
without that failure mode.)

Read-only reuse of ``bg_subtract.pipeline``'s file/ROI I/O helpers
(``_sorted_tifs``, ``_sorted_roi_files``, ``_load_rois_from_file``,
``_load_stack``, ``_save_stack``, ``_roi_mask_for_image``,
``_roi_mask_union``, ``_subtract_per_slice``, ``discover_cropped_folders``)
— does not import, modify, or depend on ``subtract_sample_folder`` (the
mean-based implementation), so the existing pipeline and its output are
completely untouched by this script.

Input layout (same as bg_subtract.py, produced by roi_cropping/)::

    <sample>/Cropped/
        gfp/gfp1.tif, gfp2.tif, ...
        cy/cy1.tif, cy2.tif, ...
        roi.zip

Output (a *different* folder name from bg_subtract.py's Background_Subtracted/,
so both methods' output can sit side by side for comparison)::

    <sample>/Background_Subtracted_Percentile/
        gfp/gfp1.tif, gfp2.tif, ...
        cy/cy1.tif, cy2.tif, ...

Two tuning knobs, both defaulting to bg_subtract.py's own defaults so the
starting point is a minimal, single-variable change:

- ``--percentile`` (default 10): which percentile of the ROI to use.
  Lower = closer to the dimmest background pixels (more conservative,
  subtracts less); higher moves back toward mean-like behavior and
  re-exposes the original problem as it approaches ~50.
- ``--gfp-multiplier`` / ``--cy-multiplier`` (default 1.25 each, matching
  bg_subtract.py): same clamp-and-scale as the existing method.

Usage::

    python bg_subtract_percentile.py --input-dir /path/to/main_folder --percentile 10
    python bg_subtract_percentile.py --input-dir /path/to/Cropped --direct --percentile 10
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import roifile

# Make the src/ layout importable without requiring `pip install -e .`.
_SRC = Path(__file__).parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from bg_subtract.pipeline import (  # noqa: E402
    _load_rois_from_file,
    _load_stack,
    _roi_mask_for_image,
    _roi_mask_union,
    _save_stack,
    _sorted_roi_files,
    _sorted_tifs,
    _subtract_per_slice,
    discover_cropped_folders,
)

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------


def compute_per_slice_subtract_values_percentile(
    stack: np.ndarray,
    mask: np.ndarray,
    percentile: float,
    multiplier: float,
    *,
    floor: float = 100.0,
    ceiling: float = 5000.0,
) -> np.ndarray:
    """Per-slice subtraction value: ``percentile(pixels inside mask) * multiplier``, clamped.

    Same shape as bg_subtract.pipeline's mean-based version — only the
    summary statistic changes (percentile instead of mean), so a bright,
    spatially-small region inside the mask (e.g. the Golgi) can't drag the
    value up the way it drags a mean up.
    """
    if not mask.any():
        return np.full(stack.shape[0], floor, dtype=np.float32)
    values = np.array(
        [np.percentile(stack[z][mask], percentile) for z in range(stack.shape[0])],
        dtype=np.float64,
    ) * multiplier
    return np.clip(values, floor, ceiling).astype(np.float32)


# ---------------------------------------------------------------------------
# Per-sample orchestration
# ---------------------------------------------------------------------------


@dataclass
class SampleOutputs:
    gfp: list[Path]
    cy: list[Path]


def subtract_sample_folder_percentile(
    cropped_root: Path,
    *,
    gfp_dirname: str = "gfp",
    cy_dirname: str = "cy",
    roi_zip_name: str = "roi.zip",
    roi_dirname: str = "roi",
    percentile: float = 10.0,
    gfp_multiplier: float = 1.25,
    cy_multiplier: float = 1.25,
    floor: float = 100.0,
    ceiling: float = 5000.0,
    output_dirname: str = "Background_Subtracted_Percentile",
) -> SampleOutputs:
    """Percentile-based background subtraction for one Cropped/-shaped folder.

    Same ROI-matching contract as bg_subtract.py's subtract_sample_folder
    (single roi.zip, one ROI per image by index; or a roi/ folder of
    per-image ROI files). Output goes to a sibling ``<output_dirname>/``
    (default ``Background_Subtracted_Percentile``, distinct from
    bg_subtract.py's ``Background_Subtracted``, so both can coexist).
    """
    cropped_root = Path(cropped_root)
    gfp_dir = cropped_root / gfp_dirname
    cy_dir = cropped_root / cy_dirname
    roi_zip = cropped_root / roi_zip_name
    roi_dir = cropped_root / roi_dirname

    if not gfp_dir.is_dir():
        raise NotADirectoryError(gfp_dir)
    if not cy_dir.is_dir():
        raise NotADirectoryError(cy_dir)

    gfp_files = _sorted_tifs(gfp_dir)
    cy_files = _sorted_tifs(cy_dir)
    if len(gfp_files) != len(cy_files):
        raise ValueError(
            f"GFP/CY count mismatch under {cropped_root}: "
            f"gfp={len(gfp_files)}, cy={len(cy_files)}"
        )
    if not gfp_files:
        log.info("%s: no images to process", cropped_root)
        return SampleOutputs([], [])

    rois: list[roifile.ImagejRoi] | None = None
    roi_files: list[Path] | None = None
    if roi_zip.is_file():
        rois = roifile.roiread(str(roi_zip))
        if isinstance(rois, roifile.ImagejRoi):
            rois = [rois]
        if len(rois) < len(gfp_files):
            raise ValueError(
                f"Not enough ROIs in {roi_zip}: rois={len(rois)}, "
                f"images={len(gfp_files)}"
            )
    elif roi_dir.is_dir():
        roi_files = _sorted_roi_files(roi_dir)
        if len(roi_files) < len(gfp_files):
            raise ValueError(
                f"Not enough ROI files in {roi_dir}: rois={len(roi_files)}, "
                f"images={len(gfp_files)}"
            )
    else:
        raise FileNotFoundError(
            f"Neither {roi_zip} nor {roi_dir}/ found — need one ROI per "
            f"image as either a single {roi_zip_name} or a {roi_dirname}/ "
            "folder of per-image ROI files."
        )

    out_root = cropped_root.parent / output_dirname
    out_gfp_dir = out_root / "gfp"
    out_cy_dir = out_root / "cy"
    out_gfp_dir.mkdir(parents=True, exist_ok=True)
    out_cy_dir.mkdir(parents=True, exist_ok=True)

    out_gfp: list[Path] = []
    out_cy: list[Path] = []
    for i, (gfp_path, cy_path) in enumerate(zip(gfp_files, cy_files)):
        gfp_stack = _load_stack(gfp_path)
        cy_stack = _load_stack(cy_path)
        hw = gfp_stack.shape[1:]
        if rois is not None:
            mask = _roi_mask_for_image(rois[i], hw)
        else:
            mask = _roi_mask_union(_load_rois_from_file(roi_files[i]), hw)  # type: ignore[index]

        gfp_sub = compute_per_slice_subtract_values_percentile(
            gfp_stack, mask, percentile, gfp_multiplier, floor=floor, ceiling=ceiling,
        )
        cy_sub = compute_per_slice_subtract_values_percentile(
            cy_stack, mask, percentile, cy_multiplier, floor=floor, ceiling=ceiling,
        )
        log.info(
            "[%s] pair %d (%s / %s): gfp_sub=[%s], cy_sub=[%s]",
            cropped_root.parent.name, i + 1, gfp_path.name, cy_path.name,
            ", ".join(f"{v:.1f}" for v in gfp_sub),
            ", ".join(f"{v:.1f}" for v in cy_sub),
        )

        gfp_out_arr = _subtract_per_slice(gfp_stack, gfp_sub)
        cy_out_arr = _subtract_per_slice(cy_stack, cy_sub)

        gfp_out_path = out_gfp_dir / gfp_path.name
        cy_out_path = out_cy_dir / cy_path.name
        _save_stack(gfp_out_path, gfp_out_arr, gfp_stack.dtype)
        _save_stack(cy_out_path, cy_out_arr, cy_stack.dtype)
        out_gfp.append(gfp_out_path)
        out_cy.append(cy_out_path)

    return SampleOutputs(out_gfp, out_cy)


def run_pipeline(
    input_dir: Path,
    *,
    roi_zip_name: str = "roi.zip",
    roi_dirname: str = "roi",
    percentile: float = 10.0,
    gfp_multiplier: float = 1.25,
    cy_multiplier: float = 1.25,
    floor: float = 100.0,
    ceiling: float = 5000.0,
    output_dirname: str = "Background_Subtracted_Percentile",
) -> dict[Path, SampleOutputs]:
    """Walk input_dir for Cropped/ folders and percentile-subtract each."""
    input_dir = Path(input_dir)
    if not input_dir.is_dir():
        raise NotADirectoryError(input_dir)

    cropped = discover_cropped_folders(
        input_dir, roi_zip_name=roi_zip_name, roi_dirname=roi_dirname,
    )
    if not cropped:
        log.warning("No Cropped/ folders found under %s", input_dir)
        return {}

    results: dict[Path, SampleOutputs] = {}
    for cropped_root in cropped:
        log.info("=== Percentile background-subtracting: %s ===", cropped_root)
        results[cropped_root] = subtract_sample_folder_percentile(
            cropped_root,
            roi_zip_name=roi_zip_name,
            roi_dirname=roi_dirname,
            percentile=percentile,
            gfp_multiplier=gfp_multiplier,
            cy_multiplier=cy_multiplier,
            floor=floor,
            ceiling=ceiling,
            output_dirname=output_dirname,
        )
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="bg_subtract_percentile",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input-dir", type=Path, required=True,
        help="Folder to walk recursively for Cropped/ subfolders. In "
        "--direct mode, this is instead the single folder to process "
        "(must directly contain gfp/, cy/, and an ROI source).",
    )
    p.add_argument(
        "--direct", action="store_true",
        help="Treat --input-dir as the sample folder itself instead of "
        "walking for Cropped/ subfolders.",
    )
    p.add_argument(
        "--roi-zip", default="roi.zip",
        help="ROI zip filename to look for inside each sample (default: roi.zip).",
    )
    p.add_argument(
        "--roi-dirname", default="roi",
        help="Fallback ROI folder name if --roi-zip isn't found (default: roi).",
    )
    p.add_argument(
        "--percentile", type=float, default=10.0,
        help="ROI percentile used in place of the mean (default: 10). "
        "Lower = more conservative (subtracts less); higher moves back "
        "toward mean-like behavior.",
    )
    p.add_argument("--gfp-multiplier", type=float, default=1.25,
                   help="Multiplier on the GFP per-cell percentile (default: 1.25).")
    p.add_argument("--cy-multiplier", type=float, default=1.25,
                   help="Multiplier on the CY per-cell percentile (default: 1.25).")
    p.add_argument("--floor", type=float, default=100.0,
                   help="Lower clamp on the subtraction value (default: 100).")
    p.add_argument("--ceiling", type=float, default=5000.0,
                   help="Upper clamp on the subtraction value (default: 5000).")
    p.add_argument(
        "--output-dirname", default="Background_Subtracted_Percentile",
        help="Output folder name, written as a sibling of Cropped/ "
        "(default: Background_Subtracted_Percentile).",
    )
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
        print(f"Percentile: {args.percentile}  multiplier(gfp/cy): "
              f"{args.gfp_multiplier}/{args.cy_multiplier}")
        outs = subtract_sample_folder_percentile(
            sample_dir,
            roi_zip_name=args.roi_zip,
            roi_dirname=args.roi_dirname,
            percentile=args.percentile,
            gfp_multiplier=args.gfp_multiplier,
            cy_multiplier=args.cy_multiplier,
            floor=args.floor,
            ceiling=args.ceiling,
            output_dirname=args.output_dirname,
        )
        out_root = sample_dir.parent / args.output_dirname
        print(f"\nDone. Background-subtracted {len(outs.gfp)} pair(s) -> {out_root}")
        return 0

    results = run_pipeline(
        args.input_dir,
        roi_zip_name=args.roi_zip,
        roi_dirname=args.roi_dirname,
        percentile=args.percentile,
        gfp_multiplier=args.gfp_multiplier,
        cy_multiplier=args.cy_multiplier,
        floor=args.floor,
        ceiling=args.ceiling,
        output_dirname=args.output_dirname,
    )
    total_pairs = sum(len(o.gfp) for o in results.values())
    print(f"\nDone. Background-subtracted {total_pairs} pair(s) "
          f"across {len(results)} sample(s).")
    for cropped_root, outs in results.items():
        out_root = cropped_root.parent / args.output_dirname
        print(f"  {cropped_root.parent}: {len(outs.gfp)} pair(s) -> {out_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
