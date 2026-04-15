#!/usr/bin/env python3
"""
Standalone colocalization analysis using NumPy.

Computes per-slice auto-thresholds (IsoData/Default), Pearson's correlation,
and thresholded Manders' coefficients for paired GFP/CY image stacks.

This is a standalone replacement for the ImageJ/JACoP workflow that runs
from the terminal.  NumPy vectorized operations make it orders of magnitude
faster than pixel-by-pixel loops in ImageJ/Jython.

Requirements:
    pip install numpy tifffile

Usage:
    python standalone_analysis.py /path/to/Background_Subtracted

Output:
    coloc_results.csv         -- per-slice Pearson's r, Manders' M1 & M2
    thresholds_multislice.csv -- per-slice thresholds (backward-compatible)
"""

import os
import sys
import csv
import argparse
import logging
from pathlib import Path
from typing import List, Optional, Union

import numpy as np

try:
    import tifffile
except ImportError:
    print("Error: tifffile is required.  Install with:  pip install tifffile")
    sys.exit(1)


log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Thresholding
# ---------------------------------------------------------------------------

def _isodata_from_histogram(hist):
    """
    Core IsoData threshold algorithm.

    Exact reimplementation of ImageJ's AutoThresholder.IJDefault():
    walks a moving index from the lowest non-zero bin, computing the mean
    of the lower and upper partitions, and stops when the index crosses the
    average of the two means.
    """
    n = len(hist)
    hist = np.asarray(hist, dtype=np.float64)

    nonzero = np.nonzero(hist)[0]
    if len(nonzero) < 2:
        return n // 2

    min_bin = int(nonzero[0])
    max_bin = int(nonzero[-1])
    if min_bin >= max_bin:
        return n // 2

    indices = np.arange(n, dtype=np.float64)
    moving = min_bin
    result = 0.0

    while True:
        # Lower partition [min_bin .. moving]
        lo = hist[min_bin:moving + 1]
        sum2 = lo.sum()
        sum1 = (indices[min_bin:moving + 1] * lo).sum()

        # Upper partition [moving+1 .. max_bin]
        hi = hist[moving + 1:max_bin + 1]
        sum4 = hi.sum()
        sum3 = (indices[moving + 1:max_bin + 1] * hi).sum()

        if sum2 == 0 or sum4 == 0:
            break

        result = (sum1 / sum2 + sum3 / sum4) / 2.0
        moving += 1

        if not ((moving + 1) <= result and moving < max_bin - 1):
            break

    return int(round(result))


def isodata_threshold(image):
    """
    Compute auto-threshold using ImageJ's Default (IsoData) method.

    Handles uint8, uint16, and float images.  Falls back to
    mean + 2*stddev when IsoData returns 0.
    """
    if image.dtype == np.uint8:
        hist = np.bincount(image.ravel(), minlength=256)
    elif image.dtype == np.uint16:
        hist = np.bincount(image.ravel(), minlength=65536)
    else:
        # Float images: quantize into 256 bins, map back afterwards
        img_min, img_max = float(image.min()), float(image.max())
        if img_min == img_max:
            return img_min
        hist, _ = np.histogram(image.ravel(), bins=256,
                               range=(img_min, img_max))
        thr_idx = _isodata_from_histogram(hist)
        return img_min + thr_idx * (img_max - img_min) / 256.0

    thr = _isodata_from_histogram(hist)
    if thr <= 0:
        return float(image.mean() + 2.0 * image.std())
    return float(thr)


# ---------------------------------------------------------------------------
# Colocalization metrics
# ---------------------------------------------------------------------------

def compute_pearson(a, b):
    """Pearson's correlation coefficient over all pixels."""
    a = a.astype(np.float64).ravel()
    b = b.astype(np.float64).ravel()

    a_dev = a - a.mean()
    b_dev = b - b.mean()

    denom = np.sqrt(np.sum(a_dev ** 2) * np.sum(b_dev ** 2))
    if denom == 0:
        return 0.0
    return float(np.sum(a_dev * b_dev) / denom)


def compute_manders(a, b, thr_a, thr_b):
    """
    Thresholded Manders' colocalization coefficients.

    M1 = sum(A | A > thrA AND B > thrB) / sum(A | A > thrA)
    M2 = sum(B | A > thrA AND B > thrB) / sum(B | B > thrB)
    """
    a = a.astype(np.float64).ravel()
    b = b.astype(np.float64).ravel()

    mask_a = a > thr_a
    mask_b = b > thr_b
    coloc = mask_a & mask_b

    sum_a_above = a[mask_a].sum()
    sum_b_above = b[mask_b].sum()

    m1 = float(a[coloc].sum() / sum_a_above) if sum_a_above > 0 else 0.0
    m2 = float(b[coloc].sum() / sum_b_above) if sum_b_above > 0 else 0.0

    return m1, m2


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------

def get_num(filename, prefix):
    """Extract number from filename after prefix."""
    base = os.path.basename(filename)
    s = base[len(prefix):-4]
    try:
        return int(s)
    except ValueError:
        return s


def get_files(folder, prefix):
    """Get sorted list of TIFF files with given prefix."""
    return sorted([f for f in os.listdir(folder)
                   if f.lower().endswith(('.tif', '.tiff'))
                   and f.lower().startswith(prefix.lower())])


def load_stack(filepath):
    """
    Load a TIFF file as a 3D array (slices, height, width).

    Single-plane images are wrapped in a leading dimension so the caller
    can always iterate over axis 0.
    """
    data = tifffile.imread(filepath)
    if data.ndim == 2:
        data = data[np.newaxis, ...]
    elif data.ndim > 3:
        # Multi-channel: keep first channel only
        if data.shape[-1] <= 4:
            data = data[..., 0]
        if data.ndim > 3:
            data = data.reshape(-1, data.shape[-2], data.shape[-1])
    return data


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Library API — reusable for batch / combined pipeline runs
# ---------------------------------------------------------------------------


def analyze_folder(
    bg_sub_dir: Union[str, Path],
    results_csv: Optional[Union[str, Path]] = None,
    thresholds_csv: Optional[Union[str, Path]] = None,
    verbose: bool = True,
) -> dict:
    """Run colocalization analysis on one ``Background_Subtracted/`` folder.

    Walks ``bg_sub_dir/gfp/`` + ``bg_sub_dir/cy/``, pairs files by the
    numeric suffix, computes per-slice IsoData thresholds, Pearson's r,
    and thresholded Manders' M1/M2, and writes two CSVs into the same
    folder (or at caller-specified paths).

    Returns a dict summary with ``pairs_processed`` and ``slices_processed``.
    """
    root = str(bg_sub_dir)
    gfp_dir = os.path.join(root, "gfp")
    cy_dir = os.path.join(root, "cy")

    if not os.path.isdir(gfp_dir):
        raise FileNotFoundError("GFP folder not found: {}".format(gfp_dir))
    if not os.path.isdir(cy_dir):
        raise FileNotFoundError("CY folder not found: {}".format(cy_dir))

    results_csv = str(results_csv) if results_csv else os.path.join(
        root, "coloc_results.csv"
    )
    thresholds_csv = str(thresholds_csv) if thresholds_csv else os.path.join(
        root, "thresholds_multislice.csv"
    )

    gfp_files = get_files(gfp_dir, "gfp")
    cy_files = get_files(cy_dir, "cy")

    say = print if verbose else log.debug
    say("GFP folder: {}".format(gfp_dir))
    say("CY folder:  {}".format(cy_dir))
    say("Results:     {}".format(results_csv))
    say("Found {} GFP and {} CY files".format(len(gfp_files), len(cy_files)))

    # ---- CSV headers ----
    with open(results_csv, 'w', newline='') as f:
        csv.writer(f).writerow(["gfp_image", "cy_image", "slice",
                                 "gfp_threshold", "cy_threshold",
                                 "pearson", "mandersA", "mandersB"])
    with open(thresholds_csv, 'w', newline='') as f:
        csv.writer(f).writerow(["gfp_image", "cy_image", "slice",
                                 "gfp_threshold", "cy_threshold"])

    total_slices = 0
    total_pairs = 0

    for gfp_file in gfp_files:
        gfp_num = get_num(gfp_file, "gfp")
        cy_file = next((f for f in cy_files if get_num(f, "cy") == gfp_num), None)

        if not cy_file:
            say("No matching CY file for {}".format(gfp_file))
            continue

        say("\nProcessing: {} + {}".format(gfp_file, cy_file))

        gfp_stack = load_stack(os.path.join(gfp_dir, gfp_file))
        cy_stack = load_stack(os.path.join(cy_dir, cy_file))

        nslices = min(gfp_stack.shape[0], cy_stack.shape[0])
        say("  {} slices".format(nslices))

        result_rows = []
        threshold_rows = []

        for s in range(nslices):
            gfp_slice = gfp_stack[s]
            cy_slice = cy_stack[s]

            thr_g = isodata_threshold(gfp_slice)
            thr_c = isodata_threshold(cy_slice)

            pearson = compute_pearson(gfp_slice, cy_slice)
            m1, m2 = compute_manders(gfp_slice, cy_slice, thr_g, thr_c)

            slice_num = s + 1
            say("  Slice {:3d}: thr_g={:7.1f}  thr_c={:7.1f}  "
                "r={:.4f}  M1={:.4f}  M2={:.4f}".format(
                    slice_num, thr_g, thr_c, pearson, m1, m2))

            result_rows.append([gfp_file, cy_file, slice_num,
                                "{:.2f}".format(thr_g), "{:.2f}".format(thr_c),
                                "{:.6f}".format(pearson),
                                "{:.6f}".format(m1), "{:.6f}".format(m2)])
            threshold_rows.append([gfp_file, cy_file, slice_num,
                                   "{:.2f}".format(thr_g), "{:.2f}".format(thr_c)])

        with open(results_csv, 'a', newline='') as f:
            csv.writer(f).writerows(result_rows)
        with open(thresholds_csv, 'a', newline='') as f:
            csv.writer(f).writerows(threshold_rows)

        total_slices += nslices
        total_pairs += 1

    say("\n" + "=" * 50)
    say("Complete!")
    say("Processed {} image pairs, {} total slices".format(total_pairs, total_slices))
    say("Results:     {}".format(results_csv))
    say("Thresholds:  {}".format(thresholds_csv))

    return {
        "bg_sub_dir": root,
        "results_csv": results_csv,
        "thresholds_csv": thresholds_csv,
        "pairs_processed": total_pairs,
        "slices_processed": total_slices,
    }


def discover_bg_sub_folders(root: Union[str, Path]) -> List[Path]:
    """Find every ``Background_Subtracted/`` folder under ``root`` that has
    both ``gfp/`` and ``cy/`` subfolders — i.e. every folder ready for
    colocalization analysis."""
    root = Path(root)
    if not root.is_dir():
        raise NotADirectoryError(root)
    found: List[Path] = []
    for path in root.rglob("Background_Subtracted"):
        if not path.is_dir():
            continue
        if not (path / "gfp").is_dir() or not (path / "cy").is_dir():
            continue
        found.append(path)
    return sorted(found)


def run_pipeline(root: Union[str, Path], verbose: bool = True) -> dict:
    """Walk ``root`` for every ``Background_Subtracted/`` folder and
    analyze each. Returns ``{folder: analysis_summary}``."""
    folders = discover_bg_sub_folders(root)
    results: dict = {}
    for folder in folders:
        if verbose:
            print("\n>>> Analyzing {}".format(folder))
        results[folder] = analyze_folder(folder, verbose=verbose)
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Colocalization analysis: Pearson's r + thresholded Manders' M1/M2"
    )
    parser.add_argument(
        "bg_sub_dir",
        help="Path to Background_Subtracted folder (must contain gfp/ and cy/ subfolders)"
    )
    args = parser.parse_args()

    try:
        analyze_folder(args.bg_sub_dir)
    except FileNotFoundError as err:
        print("Error: {}".format(err))
        sys.exit(1)


if __name__ == "__main__":
    main()
