"""Measure GFP signal inside vs. outside a Cy5-derived Golgi mask.

For each paired GFP/Cy5 cell stack the advisor wants to know how much of the
GFP signal sits within the Golgi compartment (marked by Cy5) versus outside.

Algorithm (per cell, per Z slice):
  1. Rasterise the cell ROI into a boolean *cell mask*.
  2. IsoData-threshold the Cy5 slice → binary *Golgi mask* (restricted to cell).
  3. Measure GFP pixel values inside and outside the Golgi mask, within the cell.

Input layout (any stage that has ``gfp/``, ``cy/``, and ``roi.zip``)::

    <sample>/
        gfp/gfp01.tif  gfp02.tif  …    (multi-page uint16 stacks, one per cell)
        cy/cy01.tif    cy02.tif   …    (matched 1-to-1 with gfp/)
        roi.zip                         (ImageJ ROI zip, one ROI per cell)

Works with ``Deconvoluted/``, ``Cropped/``, or ``Background_Subtracted/``
sub-folders — wherever the three ingredients coexist.

Output (written to ``<sample>/Results/``)::

    golgi_signal_analysis.csv          one row per (cell × Z slice)
    masks/<cell>_golgi_mask.tif        binary uint8 stack: 1 = Golgi, 0 = outside
    masks/<cell>_golgi_rois.zip        per-slice polygon ROIs for FIJI ROI Manager

Usage::

    python golgi_signal_analysis.py --input-dir /path/to/experiment_root

The walker finds every sub-folder under ``--input-dir`` that contains both
a ``gfp/`` (or custom) directory and a ``roi.zip`` file.
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from pathlib import Path

import numpy as np
import roifile
import tifffile
from skimage.draw import polygon as _polygon_rasterize
from skimage.measure import find_contours

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

_NUM_RE = re.compile(r"(\d+)")


def _numeric_key(name: str) -> tuple[int, str]:
    m = _NUM_RE.search(name)
    return (int(m.group(1)) if m else 10 ** 9, name.lower())


def _sorted_tifs(folder: Path) -> list[Path]:
    return sorted(
        [p for p in folder.iterdir() if p.suffix.lower() in (".tif", ".tiff")],
        key=lambda p: _numeric_key(p.name),
    )


def _load_stack(path: Path) -> np.ndarray:
    """Load a TIFF as a (Z, H, W) array, wrapping single 2-D images."""
    arr = tifffile.imread(str(path))
    if arr.ndim == 2:
        arr = arr[np.newaxis, ...]
    return arr


def _roi_polygon_xy(roi: roifile.ImagejRoi) -> np.ndarray:
    coords = roi.coordinates()
    if coords is not None and len(coords) >= 3:
        return np.asarray(coords, dtype=np.float64)
    return np.array(
        [
            (roi.left, roi.top),
            (roi.right, roi.top),
            (roi.right, roi.bottom),
            (roi.left, roi.bottom),
        ],
        dtype=np.float64,
    )


def _roi_mask(roi: roifile.ImagejRoi, image_hw: tuple[int, int]) -> np.ndarray:
    """Rasterise an ImageJ ROI polygon into a (H, W) boolean mask."""
    h, w = image_hw
    poly = _roi_polygon_xy(roi)
    rr, cc = _polygon_rasterize(poly[:, 1], poly[:, 0], shape=(h, w))
    mask = np.zeros((h, w), dtype=bool)
    mask[rr, cc] = True
    return mask


# ---------------------------------------------------------------------------
# IsoData threshold — inline port of ImageJ's Default/IsoData algorithm.
# Matches the implementation in manders_mcc/standalone_analysis.py.
# ---------------------------------------------------------------------------

def _isodata_from_histogram(hist: np.ndarray) -> int:
    hist = np.asarray(hist, dtype=np.float64)
    n = len(hist)
    nonzero = np.nonzero(hist)[0]
    if len(nonzero) < 2:
        return n // 2
    min_bin, max_bin = int(nonzero[0]), int(nonzero[-1])
    if min_bin >= max_bin:
        return n // 2
    indices = np.arange(n, dtype=np.float64)
    moving = min_bin
    result = 0.0
    while True:
        lo = hist[min_bin:moving + 1]
        hi = hist[moving + 1:max_bin + 1]
        sum_lo, sum_hi = lo.sum(), hi.sum()
        if sum_lo == 0 or sum_hi == 0:
            break
        result = (
            (indices[min_bin:moving + 1] * lo).sum() / sum_lo
            + (indices[moving + 1:max_bin + 1] * hi).sum() / sum_hi
        ) / 2.0
        moving += 1
        if not ((moving + 1) <= result and moving < max_bin - 1):
            break
    return int(round(result))


def _isodata_threshold(image: np.ndarray) -> float:
    """IsoData auto-threshold matching ImageJ's Default method."""
    if image.dtype == np.uint8:
        hist = np.bincount(image.ravel(), minlength=256)
    elif image.dtype == np.uint16:
        hist = np.bincount(image.ravel(), minlength=65536)
    else:
        img_min, img_max = float(image.min()), float(image.max())
        if img_min == img_max:
            return img_min
        hist, _ = np.histogram(image.ravel(), bins=256, range=(img_min, img_max))
        idx = _isodata_from_histogram(hist)
        return img_min + idx * (img_max - img_min) / 256.0
    thr = _isodata_from_histogram(hist)
    if thr <= 0:
        return float(image.mean() + 2.0 * image.std())
    return float(thr)


# ---------------------------------------------------------------------------
# Core measurement
# ---------------------------------------------------------------------------

CSV_COLUMNS = [
    "cell",
    "slice",
    "cy5_threshold",
    "area_cell",
    "area_inside",
    "area_outside",
    "gfp_mean_inside",
    "gfp_mean_outside",
    "gfp_integrated_inside",
    "gfp_integrated_outside",
    "gfp_integrated_total",
    "fraction_gfp_inside",
]


def measure_cell(
    gfp_stack: np.ndarray,
    cy_stack: np.ndarray,
    cell_roi: roifile.ImagejRoi,
    cell_name: str,
) -> tuple[list[dict], np.ndarray]:
    """Measure per-slice GFP signal inside vs. outside the Cy5 Golgi mask.

    Parameters
    ----------
    gfp_stack:
        (Z, H, W) array of GFP signal.
    cy_stack:
        (Z, H, W) array of Cy5 signal (Golgi marker).
    cell_roi:
        ImageJ ROI defining the cell boundary.
    cell_name:
        Identifier string (used for the ``cell`` column in the CSV).

    Returns
    -------
    rows:
        List of dicts, one per Z slice, ready to write to CSV.
    mask_stack:
        (Z, H, W) uint8 binary array — 1 inside Golgi mask, 0 outside.
    """
    hw = gfp_stack.shape[1:]
    cell_mask = _roi_mask(cell_roi, hw)  # (H, W) bool — cell boundary

    n_slices = min(gfp_stack.shape[0], cy_stack.shape[0])
    rows: list[dict] = []
    mask_stack = np.zeros((n_slices, *hw), dtype=np.uint8)

    for z in range(n_slices):
        cy_slice = cy_stack[z]
        gfp_slice = gfp_stack[z]

        threshold = _isodata_threshold(cy_slice[cell_mask])  # threshold within cell only
        golgi_mask = (cy_slice >= threshold) & cell_mask
        outside_mask = cell_mask & ~golgi_mask

        mask_stack[z] = golgi_mask.astype(np.uint8)

        area_cell = int(cell_mask.sum())
        area_inside = int(golgi_mask.sum())
        area_outside = int(outside_mask.sum())

        gfp_f = gfp_slice.astype(np.float64)

        if area_inside > 0:
            gfp_inside = gfp_f[golgi_mask]
            mean_inside = float(gfp_inside.mean())
            intg_inside = float(gfp_inside.sum())
        else:
            mean_inside = 0.0
            intg_inside = 0.0

        if area_outside > 0:
            gfp_outside = gfp_f[outside_mask]
            mean_outside = float(gfp_outside.mean())
            intg_outside = float(gfp_outside.sum())
        else:
            mean_outside = 0.0
            intg_outside = 0.0

        intg_total = intg_inside + intg_outside
        fraction_inside = intg_inside / intg_total if intg_total > 0 else 0.0

        rows.append({
            "cell": cell_name,
            "slice": z + 1,
            "cy5_threshold": f"{threshold:.2f}",
            "area_cell": area_cell,
            "area_inside": area_inside,
            "area_outside": area_outside,
            "gfp_mean_inside": f"{mean_inside:.4f}",
            "gfp_mean_outside": f"{mean_outside:.4f}",
            "gfp_integrated_inside": f"{intg_inside:.2f}",
            "gfp_integrated_outside": f"{intg_outside:.2f}",
            "gfp_integrated_total": f"{intg_total:.2f}",
            "fraction_gfp_inside": f"{fraction_inside:.6f}",
        })

    return rows, mask_stack


# ---------------------------------------------------------------------------
# Mask → ImageJ ROI conversion
# ---------------------------------------------------------------------------

def _mask_to_rois(mask_stack: np.ndarray, cell_name: str) -> list[roifile.ImagejRoi]:
    """Convert a (Z, H, W) binary mask into a list of per-slice polygon ROIs.

    Each contour in each slice becomes one ImageJ polygon ROI with its
    ``position`` set to the 1-based slice index so FIJI's ROI Manager places
    it on the correct Z plane.
    """
    rois: list[roifile.ImagejRoi] = []
    n_slices = mask_stack.shape[0]
    for z in range(n_slices):
        slice_mask = mask_stack[z].astype(bool)
        contours = find_contours(slice_mask, level=0.5)
        for c_idx, contour in enumerate(contours):
            # find_contours returns (row, col) = (y, x); roifile wants (x, y)
            xy = contour[:, ::-1].astype(np.float32)  # (N, 2) as (x, y)
            name = (
                f"{cell_name}_z{z + 1:02d}"
                if len(contours) == 1
                else f"{cell_name}_z{z + 1:02d}_r{c_idx + 1}"
            )
            roi = roifile.ImagejRoi.frompoints(
                xy,
                name=name,
                position=z + 1,  # 1-based Z slice for FIJI ROI Manager
            )
            rois.append(roi)
    return rois


# ---------------------------------------------------------------------------
# Sample-level processing
# ---------------------------------------------------------------------------

def process_sample(
    sample_dir: Path,
    *,
    gfp_dirname: str = "gfp",
    cy_dirname: str = "cy",
    roi_zip_name: str = "roi.zip",
) -> Path | None:
    """Run Golgi signal analysis for one sample folder.

    Returns the path to the output CSV, or ``None`` if nothing was processed.
    """
    roi_zip = sample_dir / roi_zip_name
    if not roi_zip.is_file():
        log.warning("No %s in %s — skipping", roi_zip_name, sample_dir)
        return None

    gfp_dir = sample_dir / gfp_dirname
    cy_dir = sample_dir / cy_dirname
    if not gfp_dir.is_dir():
        log.warning("No %s/ folder in %s — skipping", gfp_dirname, sample_dir)
        return None
    if not cy_dir.is_dir():
        log.warning("No %s/ folder in %s — skipping", cy_dirname, sample_dir)
        return None

    rois = roifile.roiread(str(roi_zip))
    if isinstance(rois, roifile.ImagejRoi):
        rois = [rois]
    if not rois:
        log.warning("Empty ROI zip %s — skipping", roi_zip)
        return None

    gfp_files = _sorted_tifs(gfp_dir)
    cy_files = _sorted_tifs(cy_dir)

    n = min(len(gfp_files), len(cy_files), len(rois))
    if n == 0:
        log.warning("%s: no matched GFP/Cy/ROI triplets — skipping", sample_dir.name)
        return None
    if len(gfp_files) != len(cy_files) or len(gfp_files) != len(rois):
        log.warning(
            "%s: GFP=%d, Cy=%d, ROIs=%d — processing first %d",
            sample_dir.name, len(gfp_files), len(cy_files), len(rois), n,
        )

    out_dir = sample_dir / "Results"
    masks_dir = out_dir / "masks"
    out_dir.mkdir(parents=True, exist_ok=True)
    masks_dir.mkdir(exist_ok=True)

    all_rows: list[dict] = []

    for i in range(n):
        gfp_stack = _load_stack(gfp_files[i])
        cy_stack = _load_stack(cy_files[i])
        cell_name = gfp_files[i].stem  # e.g. "gfp01"

        rows, mask_stack = measure_cell(gfp_stack, cy_stack, rois[i], cell_name)
        all_rows.extend(rows)

        # Binary mask TIFF
        mask_path = masks_dir / f"{cell_name}_golgi_mask.tif"
        tifffile.imwrite(
            str(mask_path),
            mask_stack,
            photometric="minisblack",
            compression=None,
        )

        # ImageJ ROI zip
        cell_rois = _mask_to_rois(mask_stack, cell_name)
        if cell_rois:
            roi_out = masks_dir / f"{cell_name}_golgi_rois.zip"
            roifile.roiwrite(str(roi_out), cell_rois)

        log.info(
            "  %s: %d slice(s), mask written to %s",
            cell_name, len(rows), masks_dir.name,
        )

    if not all_rows:
        return None

    csv_path = out_dir / "golgi_signal_analysis.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(all_rows)

    log.info("Wrote %d rows → %s", len(all_rows), csv_path)
    return csv_path


# ---------------------------------------------------------------------------
# Auto-discovery
# ---------------------------------------------------------------------------

def discover_sample_folders(
    root: Path,
    gfp_dirname: str,
    cy_dirname: str,
    roi_zip_name: str,
) -> list[Path]:
    """Find folders containing gfp/, cy/, and roi.zip anywhere under root."""
    root = Path(root)
    samples: list[Path] = []
    seen: set[Path] = set()
    for path in root.rglob(roi_zip_name):
        if not path.is_file():
            continue
        folder = path.parent
        if folder in seen:
            continue
        has_gfp = (folder / gfp_dirname).is_dir()
        has_cy = (folder / cy_dirname).is_dir()
        if has_gfp and has_cy:
            # Skip previously written output folders
            if "Results" in folder.parts:
                continue
            samples.append(folder)
            seen.add(folder)
    return sorted(samples)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Root folder to search for sample directories (recursively).",
    )
    parser.add_argument(
        "--gfp-dirname", default="gfp",
        help="GFP channel folder name (default: gfp).",
    )
    parser.add_argument(
        "--cy-dirname", default="cy",
        help="Cy5 channel folder name (default: cy).",
    )
    parser.add_argument(
        "--roi-zip", default="roi.zip",
        help="ROI zip filename to look for inside each sample (default: roi.zip).",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose logging.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stdout,
    )

    root = args.input_dir.resolve()
    if not root.is_dir():
        parser.error(f"--input-dir is not a directory: {root}")

    samples = discover_sample_folders(
        root, args.gfp_dirname, args.cy_dirname, args.roi_zip,
    )

    if not samples:
        print(
            f"No sample folders found under {root}\n"
            f"(looking for folders with {args.gfp_dirname}/, "
            f"{args.cy_dirname}/, and {args.roi_zip})"
        )
        sys.exit(1)

    print(f"Found {len(samples)} sample folder(s)")
    total = 0
    for sample in samples:
        print(f"  Processing: {sample.name}")
        csv_path = process_sample(
            sample,
            gfp_dirname=args.gfp_dirname,
            cy_dirname=args.cy_dirname,
            roi_zip_name=args.roi_zip,
        )
        if csv_path:
            total += 1

    print(f"\nDone. Wrote {total} CSV file(s).")


if __name__ == "__main__":
    main()
