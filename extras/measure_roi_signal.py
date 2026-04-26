"""Measure per-slice signal intensity inside each cell's ROI.

Python port of ``Measure_ROI_Signal_Per_Slice.ijm``. For every cropped
cell stack, measures the mean, integrated density, area, min, max, and
standard deviation of pixels inside the ROI on each Z slice, then writes
a single combined CSV with all cells and slices.

The original ImageJ macro saved one CSV per cell. This script writes one
CSV per sample (all cells, all slices) which is easier to work with in
downstream analysis (R, pandas, etc.).

Input layout (output of ``roi_cropping/`` or ``background_subtraction/``)::

    <sample>/
        Cropped/           (or Background_Subtracted/)
            gfp/gfp01.tif
            cy/cy01.tif
            roi.zip

Output::

    <sample>/
        Results/
            signal_measurements.csv

Usage::

    python measure_roi_signal.py --input-dir /path/to/main_folder

The walker finds every folder containing ``{gfp,cy}/`` + ``roi.zip``
anywhere under ``--input-dir``.
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import roifile
import tifffile
from skimage.draw import polygon as _polygon_rasterize

log = logging.getLogger(__name__)

_NUM_RE = re.compile(r"(\d+)")


def _numeric_key(name: str) -> tuple[int, str]:
    m = _NUM_RE.search(name)
    return (int(m.group(1)) if m else 10**9, name.lower())


def _sorted_tifs(folder: Path) -> list[Path]:
    return sorted(
        [p for p in folder.iterdir() if p.suffix.lower() in (".tif", ".tiff")],
        key=lambda p: _numeric_key(p.name),
    )


def _load_stack(path: Path) -> np.ndarray:
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
    h, w = image_hw
    poly = _roi_polygon_xy(roi)
    rr, cc = _polygon_rasterize(poly[:, 1], poly[:, 0], shape=(h, w))
    mask = np.zeros((h, w), dtype=bool)
    mask[rr, cc] = True
    return mask


@dataclass
class SliceMeasurement:
    cell: str
    channel: str
    image_file: str
    slice_index: int
    area: int
    mean: float
    std: float
    min_val: float
    max_val: float
    integrated_density: float


def measure_stack(
    stack: np.ndarray,
    roi: roifile.ImagejRoi,
    cell_name: str,
    channel: str,
    image_file: str,
) -> list[SliceMeasurement]:
    """Measure per-slice statistics inside the ROI for one stack."""
    mask = _roi_mask(roi, stack.shape[1:])
    area = int(mask.sum())
    if area == 0:
        return []

    results = []
    for z in range(stack.shape[0]):
        pixels = stack[z][mask].astype(np.float64)
        results.append(SliceMeasurement(
            cell=cell_name,
            channel=channel,
            image_file=image_file,
            slice_index=z + 1,
            area=area,
            mean=float(pixels.mean()),
            std=float(pixels.std()),
            min_val=float(pixels.min()),
            max_val=float(pixels.max()),
            integrated_density=float(pixels.sum()),
        ))
    return results


CSV_COLUMNS = [
    "cell", "channel", "image_file", "slice",
    "area", "mean", "std_dev", "min", "max", "integrated_density",
]


def process_sample(
    sample_dir: Path,
    *,
    gfp_dirname: str = "gfp",
    cy_dirname: str = "cy",
    roi_zip_name: str = "roi.zip",
    channels: list[str] | None = None,
) -> Path | None:
    """Measure ROI signal for one sample folder.

    Returns the path to the output CSV, or None if nothing was processed.
    """
    if channels is None:
        channels = [gfp_dirname, cy_dirname]

    roi_zip = sample_dir / roi_zip_name
    if not roi_zip.is_file():
        log.warning("No %s in %s — skipping", roi_zip_name, sample_dir)
        return None

    rois = roifile.roiread(str(roi_zip))
    if isinstance(rois, roifile.ImagejRoi):
        rois = [rois]
    if not rois:
        log.warning("Empty ROI zip %s — skipping", roi_zip)
        return None

    all_measurements: list[SliceMeasurement] = []

    for channel in channels:
        ch_dir = sample_dir / channel
        if not ch_dir.is_dir():
            log.info("Channel folder %s not found — skipping channel", ch_dir)
            continue

        tif_files = _sorted_tifs(ch_dir)
        n = min(len(tif_files), len(rois))
        if len(tif_files) != len(rois):
            log.warning(
                "%s: %s has %d TIFFs but %d ROIs — processing first %d",
                sample_dir.name, channel, len(tif_files), len(rois), n,
            )

        for i in range(n):
            stack = _load_stack(tif_files[i])
            cell_name = tif_files[i].stem
            measurements = measure_stack(
                stack, rois[i], cell_name, channel, tif_files[i].name,
            )
            all_measurements.extend(measurements)
            log.debug(
                "  %s/%s: %d slices, ROI area=%d",
                channel, tif_files[i].name, stack.shape[0],
                measurements[0].area if measurements else 0,
            )

    if not all_measurements:
        return None

    out_dir = sample_dir / "Results"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "signal_measurements.csv"

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(CSV_COLUMNS)
        for m in all_measurements:
            writer.writerow([
                m.cell, m.channel, m.image_file, m.slice_index,
                m.area, f"{m.mean:.4f}", f"{m.std:.4f}",
                f"{m.min_val:.1f}", f"{m.max_val:.1f}",
                f"{m.integrated_density:.2f}",
            ])

    log.info("Wrote %d rows to %s", len(all_measurements), csv_path)
    return csv_path


def discover_sample_folders(
    root: Path,
    gfp_dirname: str,
    cy_dirname: str,
    roi_zip_name: str,
) -> list[Path]:
    """Find folders containing channel dirs + roi.zip anywhere under root."""
    root = Path(root)
    samples: list[Path] = []
    seen: set[Path] = set()
    for path in root.rglob(roi_zip_name):
        if not path.is_file():
            continue
        folder = path.parent
        if folder in seen:
            continue
        has_channel = (folder / gfp_dirname).is_dir() or (folder / cy_dirname).is_dir()
        if has_channel:
            samples.append(folder)
            seen.add(folder)
    return sorted(samples)


def _measure_direct(
    channel_dir: Path,
    roi_zip_path: Path,
    channel_name: str,
) -> Path | None:
    """Direct mode: measure one channel folder against a specific ROI zip.

    Used when the user points at a single channel folder + ROI zip
    rather than a sample root.
    """
    if not channel_dir.is_dir():
        raise NotADirectoryError(channel_dir)
    if not roi_zip_path.is_file():
        raise FileNotFoundError(roi_zip_path)

    rois = roifile.roiread(str(roi_zip_path))
    if isinstance(rois, roifile.ImagejRoi):
        rois = [rois]
    if not rois:
        log.warning("Empty ROI zip %s", roi_zip_path)
        return None

    tif_files = _sorted_tifs(channel_dir)
    n = min(len(tif_files), len(rois))
    if not n:
        log.warning("No TIFFs in %s or no ROIs — nothing to measure", channel_dir)
        return None
    if len(tif_files) != len(rois):
        log.warning(
            "%d TIFFs but %d ROIs — processing first %d",
            len(tif_files), len(rois), n,
        )

    all_measurements: list[SliceMeasurement] = []
    for i in range(n):
        stack = _load_stack(tif_files[i])
        measurements = measure_stack(
            stack, rois[i], tif_files[i].stem, channel_name, tif_files[i].name,
        )
        all_measurements.extend(measurements)

    if not all_measurements:
        return None

    out_dir = channel_dir.parent / "Results"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "signal_measurements.csv"

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(CSV_COLUMNS)
        for m in all_measurements:
            writer.writerow([
                m.cell, m.channel, m.image_file, m.slice_index,
                m.area, f"{m.mean:.4f}", f"{m.std:.4f}",
                f"{m.min_val:.1f}", f"{m.max_val:.1f}",
                f"{m.integrated_density:.2f}",
            ])

    log.info("Wrote %d rows to %s", len(all_measurements), csv_path)
    return csv_path


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="measure_roi_signal",
        description="Measure per-slice signal intensity inside each cell's "
        "ROI. Writes a single combined CSV per sample.\n\n"
        "Two modes:\n"
        "  Auto-discovery: --input-dir points at a root folder; the script\n"
        "    finds every subfolder with gfp/cy dirs + roi.zip.\n"
        "  Direct: --input-dir points at a single channel folder (e.g.\n"
        "    .../Background_Subtracted/cy) and --roi-zip is the full path\n"
        "    to the ROI zip.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Root folder to search for sample directories, OR a single "
        "channel folder (e.g. .../Cropped/cy) when using --roi-zip "
        "with a full path.",
    )
    p.add_argument(
        "--gfp-dirname", default="gfp",
        help="GFP channel folder name (default: gfp).",
    )
    p.add_argument(
        "--cy-dirname", default="cy",
        help="Cy channel folder name (default: cy).",
    )
    p.add_argument(
        "--roi-zip", default="roi.zip",
        help="ROI zip filename to look for inside each sample (default: "
        "roi.zip). Can also be a full path to a specific zip file, in "
        "which case --input-dir is treated as a single channel folder.",
    )
    p.add_argument(
        "--channels",
        nargs="+",
        default=None,
        help="Channel folder names to measure (default: both gfp and cy). "
        "Pass just one to measure a single channel.",
    )
    p.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose logging.",
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

    root = args.input_dir.resolve()
    if not root.is_dir():
        parser.error(f"--input-dir is not a directory: {root}")

    roi_zip_arg = Path(args.roi_zip)

    # Direct mode: --roi-zip is a full path to a specific zip file.
    if roi_zip_arg.is_absolute() or roi_zip_arg.is_file():
        roi_zip_path = roi_zip_arg.resolve()
        if not roi_zip_path.is_file():
            parser.error(f"ROI zip not found: {roi_zip_path}")

        # Infer channel name from the folder name (e.g. "cy", "gfp").
        channel_name = root.name
        print(f"Direct mode: measuring {root}")
        print(f"  ROI zip: {roi_zip_path}")
        print(f"  Channel: {channel_name}")

        csv_path = _measure_direct(root, roi_zip_path, channel_name)
        if csv_path:
            print(f"\nDone. Wrote {csv_path}")
        else:
            print("\nNo measurements produced.")
        return

    # Auto-discovery mode: --roi-zip is just a filename (e.g. "roi.zip").
    roi_zip_name = args.roi_zip
    samples = discover_sample_folders(
        root, args.gfp_dirname, args.cy_dirname, roi_zip_name,
    )

    if not samples:
        print(f"No sample folders found under {root}")
        print(f"(looking for folders with {args.gfp_dirname}/ or "
              f"{args.cy_dirname}/ plus {roi_zip_name})")
        print()
        print("Tip: you can also point --input-dir at a single channel "
              "folder and pass --roi-zip as a full path:")
        print(f"  python measure_roi_signal.py --input-dir /path/to/cy "
              f"--roi-zip /path/to/roi.zip")
        sys.exit(1)

    print(f"Found {len(samples)} sample folder(s)")

    channels = args.channels or [args.gfp_dirname, args.cy_dirname]
    total_csvs = 0

    for sample in samples:
        print(f"  Processing: {sample.name}")
        csv_path = process_sample(
            sample,
            gfp_dirname=args.gfp_dirname,
            cy_dirname=args.cy_dirname,
            roi_zip_name=roi_zip_name,
            channels=channels,
        )
        if csv_path:
            total_csvs += 1

    print(f"\nDone. Wrote {total_csvs} CSV file(s).")


if __name__ == "__main__":
    main()
