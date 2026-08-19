"""Erase the ROI-enclosed region (e.g. a cell body) from paired GFP/Cy stacks.

For each GFP/Cy image pair, opens the corresponding ImageJ ROI and zeroes
out every pixel inside that ROI, on every Z slice of both channels. Useful
for isolating signal *outside* a drawn region (e.g. removing the cell body
to look at background / neurite signal).

Input layout::

    <parent_dir>/
        gfp/    gfp01.tif  gfp02.tif  ...
        cy/     cy01.tif   cy02.tif   ...
        roi/    01.zip     02.zip     ...

Files in each folder are paired positionally after sorting by the numeric
part of their filename — this matches the numbering already used by this
project's ROI-drawing output (e.g. ``roi_original/01.zip`` aligned 1:1 with
sequentially-renumbered channel stacks).

Output (written to ``<parent_dir>/Body_Removed/``)::

    Body_Removed/
        gfp/    gfp01.tif  gfp02.tif  ...
        cy/     cy01.tif   cy02.tif   ...

Usage::

    python remove_body.py --input-dir /path/to/parent_dir
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path

import numpy as np
import roifile
import tifffile
from skimage.draw import polygon as _polygon_rasterize

log = logging.getLogger(__name__)

_NUM_RE = re.compile(r"(\d+)")


def _numeric_key(name: str) -> tuple[int, str]:
    m = _NUM_RE.search(name)
    return (int(m.group(1)) if m else 10 ** 9, name.lower())


def _sorted_tifs(folder: Path) -> list[Path]:
    return sorted(
        [p for p in folder.iterdir() if p.suffix.lower() in (".tif", ".tiff")],
        key=lambda p: _numeric_key(p.name),
    )


def _sorted_roi_zips(folder: Path) -> list[Path]:
    return sorted(
        [p for p in folder.iterdir() if p.suffix.lower() == ".zip"],
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
    """Rasterise a single ImageJ ROI polygon into a (H, W) boolean mask."""
    h, w = image_hw
    poly = _roi_polygon_xy(roi)
    rr, cc = _polygon_rasterize(poly[:, 1], poly[:, 0], shape=(h, w))
    mask = np.zeros((h, w), dtype=bool)
    mask[rr, cc] = True
    return mask


def _combined_roi_mask(
    rois: list[roifile.ImagejRoi], image_hw: tuple[int, int]
) -> np.ndarray:
    """Union of every ROI in a zip into one (H, W) boolean mask.

    A ROI zip can hold more than one polygon (e.g. a cell body split into
    parts); every one of them is treated as "body" and gets erased.
    """
    mask = np.zeros(image_hw, dtype=bool)
    for roi in rois:
        mask |= _roi_mask(roi, image_hw)
    return mask


def erase_roi_from_stack(stack: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Zero out ``mask``'s True region on every Z slice of ``stack``."""
    out = stack.copy()
    out[:, mask] = 0
    return out


def _load_rois(roi_zip: Path) -> list[roifile.ImagejRoi]:
    rois = roifile.roiread(str(roi_zip))
    if isinstance(rois, roifile.ImagejRoi):
        rois = [rois]
    return list(rois)


def process_parent_folder(
    parent_dir: Path,
    *,
    gfp_dirname: str = "gfp",
    cy_dirname: str = "cy",
    roi_dirname: str = "roi",
    output_dirname: str = "Body_Removed",
) -> int:
    """Erase each image's ROI region from its GFP/Cy stacks.

    Returns the number of image pairs processed.
    """
    gfp_dir = parent_dir / gfp_dirname
    cy_dir = parent_dir / cy_dirname
    roi_dir = parent_dir / roi_dirname
    for label, d in (("GFP", gfp_dir), ("Cy", cy_dir), ("ROI", roi_dir)):
        if not d.is_dir():
            raise NotADirectoryError(f"{label} folder not found: {d}")

    gfp_files = _sorted_tifs(gfp_dir)
    cy_files = _sorted_tifs(cy_dir)
    roi_files = _sorted_roi_zips(roi_dir)

    n = min(len(gfp_files), len(cy_files), len(roi_files))
    if n == 0:
        log.warning(
            "%s: no matched GFP/Cy/ROI triplets (GFP=%d, Cy=%d, ROI=%d)",
            parent_dir, len(gfp_files), len(cy_files), len(roi_files),
        )
        return 0
    if len(gfp_files) != len(cy_files) or len(gfp_files) != len(roi_files):
        log.warning(
            "%s: GFP=%d, Cy=%d, ROI=%d — processing first %d (matched by "
            "sorted position, not filename)",
            parent_dir, len(gfp_files), len(cy_files), len(roi_files), n,
        )

    out_gfp_dir = parent_dir / output_dirname / gfp_dirname
    out_cy_dir = parent_dir / output_dirname / cy_dirname
    out_gfp_dir.mkdir(parents=True, exist_ok=True)
    out_cy_dir.mkdir(parents=True, exist_ok=True)

    for i in range(n):
        gfp_path, cy_path, roi_path = gfp_files[i], cy_files[i], roi_files[i]

        gfp_stack = _load_stack(gfp_path)
        cy_stack = _load_stack(cy_path)
        if gfp_stack.shape[1:] != cy_stack.shape[1:]:
            log.warning(
                "%s / %s: GFP shape %s != Cy shape %s — skipping",
                gfp_path.name, cy_path.name, gfp_stack.shape, cy_stack.shape,
            )
            continue

        rois = _load_rois(roi_path)
        if not rois:
            log.warning("%s: empty ROI zip — skipping", roi_path.name)
            continue

        mask = _combined_roi_mask(rois, gfp_stack.shape[1:])

        gfp_out = erase_roi_from_stack(gfp_stack, mask)
        cy_out = erase_roi_from_stack(cy_stack, mask)

        tifffile.imwrite(
            str(out_gfp_dir / gfp_path.name), gfp_out,
            imagej=True, photometric="minisblack", compression=None,
        )
        tifffile.imwrite(
            str(out_cy_dir / cy_path.name), cy_out,
            imagej=True, photometric="minisblack", compression=None,
        )

        log.info(
            "%s / %s <- %s: erased %d px region across %d slice(s)",
            gfp_path.name, cy_path.name, roi_path.name,
            int(mask.sum()), gfp_stack.shape[0],
        )

    return n


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Parent folder containing gfp/, cy/, and roi/ subfolders.",
    )
    parser.add_argument(
        "--gfp-dirname", default="gfp",
        help="GFP channel folder name (default: gfp).",
    )
    parser.add_argument(
        "--cy-dirname", default="cy",
        help="Cy channel folder name (default: cy).",
    )
    parser.add_argument(
        "--roi-dirname", default="roi",
        help="ROI folder name, containing one .zip per image (default: roi).",
    )
    parser.add_argument(
        "--output-dirname", default="Body_Removed",
        help="Output folder name, written under --input-dir (default: Body_Removed).",
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

    parent_dir = args.input_dir.resolve()
    if not parent_dir.is_dir():
        parser.error(f"--input-dir is not a directory: {parent_dir}")

    print(f"\n=== Remove body ROI ===")
    print(f"Input folder:  {parent_dir}")
    print(f"Output folder: {parent_dir / args.output_dirname}\n")

    try:
        n = process_parent_folder(
            parent_dir,
            gfp_dirname=args.gfp_dirname,
            cy_dirname=args.cy_dirname,
            roi_dirname=args.roi_dirname,
            output_dirname=args.output_dirname,
        )
    except NotADirectoryError as e:
        parser.error(str(e))
        return

    print(f"\n=== Done: processed {n} image pair(s) ===")


if __name__ == "__main__":
    main()
