"""Core single-pass cropping logic.

Replaces the MATLAB ``H_roi_duplicate_image_all_channels_Recursive.m`` +
ImageJ ``ROI_Select_Duplicate_TIFF_Loop.ijm`` two-step flow with a single
pass that emits one cropped single-cell stack per ROI directly.

Per-sample input layout (produced by ``preprocessing/`` + the manual ImageJ
``Color_Merge_Automated_PreloadROIs_Adjust`` macro)::

    <sample>/
        gfp/gfpN_decon.tif           ← multi-page stacks, one per FOV
        cy/cyN_decon.tif             ← matched 1:1 to gfp by sort order
        roi/NN.zip                   ← user-edited ROIs, one zip per FOV

Output layout::

    <sample>/Cropped/
        gfp/gfp1.tif, gfp2.tif, ...  ← one stack per ROI (sequential)
        cy/cy1.tif, cy2.tif, ...     ← matched 1:1
        roi/1.roi, 2.roi, ...        ← per-cell ROI in bbox-local coords
        roi.zip                      ← combined zip of all per-cell ROIs

Cropping reproduces the ImageJ ``Duplicate (with ROI)`` + ``Make Inverse``
+ ``Clear stack`` behavior: the output is the tight bounding box of the
ROI with all pixels outside the ROI polygon zeroed.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import roifile
import tifffile
from skimage.draw import polygon as _polygon_rasterize

log = logging.getLogger(__name__)


_NUM_RE = re.compile(r"(\d+)")


def _numeric_key(name: str) -> tuple[int, str]:
    """Sort key matching the ImageJ macro: first contiguous digit run, then name."""
    m = _NUM_RE.search(name)
    return (int(m.group(1)) if m else 10**9, name.lower())


def _sorted_files(folder: Path, pattern: str) -> list[Path]:
    return sorted(folder.glob(pattern), key=lambda p: _numeric_key(p.name))


def _load_stack(path: Path) -> np.ndarray:
    """Load a TIFF as a (Z, H, W) array. 2D files become (1, H, W)."""
    arr = tifffile.imread(str(path))
    if arr.ndim == 2:
        arr = arr[np.newaxis, ...]
    return arr


def _save_stack(path: Path, stack: np.ndarray, dtype: np.dtype) -> None:
    """Save a (Z, H, W) array as an uncompressed multi-page TIFF."""
    out = stack.astype(dtype, copy=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(str(path), out, photometric="minisblack", compression=None)


# ---------------------------------------------------------------------------
# ROI handling
# ---------------------------------------------------------------------------


def _roi_polygon_xy(roi: roifile.ImagejRoi) -> np.ndarray:
    """Return the ROI polygon as an (N, 2) array of (x, y) image coords.

    Falls back to the bbox rectangle if the ROI has no explicit polygon
    coordinates (e.g. plain rectangle ROIs).
    """
    coords = roi.coordinates()
    if coords is not None and len(coords) >= 3:
        return np.asarray(coords, dtype=np.float64)
    # Rectangle fallback — corners of (left, top) -> (right, bottom).
    return np.array(
        [
            (roi.left, roi.top),
            (roi.right, roi.top),
            (roi.right, roi.bottom),
            (roi.left, roi.bottom),
        ],
        dtype=np.float64,
    )


def _roi_bbox(roi: roifile.ImagejRoi, image_shape: tuple[int, int]) -> tuple[int, int, int, int]:
    """Return the integer bbox (x0, y0, x1, y1) clipped to the image.

    Uses the ImageJ-convention bbox stored on the ROI itself
    (``right = max(x) + 1`` etc.; half-open on the right/bottom). Falls
    back to recomputing from polygon vertices if the ROI lacks an
    explicit bbox.
    """
    h, w = image_shape
    x0 = int(roi.left)
    y0 = int(roi.top)
    x1 = int(roi.right)
    y1 = int(roi.bottom)
    if x1 <= x0 or y1 <= y0:
        # Fallback: recompute from the polygon (matches ImageJ width
        # convention of max - min + 1).
        poly = _roi_polygon_xy(roi)
        x0 = int(np.floor(poly[:, 0].min()))
        y0 = int(np.floor(poly[:, 1].min()))
        x1 = int(np.floor(poly[:, 0].max())) + 1
        y1 = int(np.floor(poly[:, 1].max())) + 1
    # Clip to image extents.
    x0 = max(0, min(w, x0))
    y0 = max(0, min(h, y0))
    x1 = max(0, min(w, x1))
    y1 = max(0, min(h, y1))
    if x1 <= x0 or y1 <= y0:
        raise ValueError(
            f"ROI '{roi.name}' has empty bbox after clipping to image "
            f"shape {image_shape}: x=[{x0},{x1}), y=[{y0},{y1})"
        )
    return x0, y0, x1, y1


def _local_polygon_mask(
    roi: roifile.ImagejRoi,
    bbox: tuple[int, int, int, int],
) -> tuple[np.ndarray, np.ndarray]:
    """Rasterize the ROI polygon into a bbox-local boolean mask.

    Returns ``(mask, local_polygon_xy)`` where ``mask`` has shape
    ``(bbox_h, bbox_w)`` and ``local_polygon_xy`` is the polygon coords
    re-origined to the bbox.
    """
    x0, y0, x1, y1 = bbox
    h_local = y1 - y0
    w_local = x1 - x0

    poly = _roi_polygon_xy(roi)
    local = poly.copy()
    local[:, 0] -= x0
    local[:, 1] -= y0

    # skimage.draw.polygon takes (rows, cols).
    rr, cc = _polygon_rasterize(local[:, 1], local[:, 0], shape=(h_local, w_local))
    mask = np.zeros((h_local, w_local), dtype=bool)
    mask[rr, cc] = True
    return mask, local


def _adjust_roi_to_bbox(
    roi: roifile.ImagejRoi,
    bbox: tuple[int, int, int, int],
) -> roifile.ImagejRoi:
    """Return a new ROI with all coords re-origined to the bbox.

    The output ROI is sized exactly to the bbox so it overlays correctly
    on the cropped image.
    """
    x0, y0, x1, y1 = bbox
    poly = _roi_polygon_xy(roi)
    local = poly.copy()
    local[:, 0] -= x0
    local[:, 1] -= y0
    new = roifile.ImagejRoi.frompoints(local)
    new.name = roi.name
    return new


# ---------------------------------------------------------------------------
# Per-pair / per-sample orchestration
# ---------------------------------------------------------------------------


@dataclass
class PairOutputs:
    gfp: list[Path]
    cy: list[Path]
    roi: list[Path]


def crop_pair_with_rois(
    gfp_path: Path,
    cy_path: Path,
    roi_zip_path: Path,
    out_gfp_dir: Path,
    out_cy_dir: Path,
    out_roi_dir: Path,
    start_index: int,
    *,
    gfp_prefix: str = "gfp",
    cy_prefix: str = "cy",
) -> tuple[PairOutputs, list[roifile.ImagejRoi]]:
    """Crop one (gfp, cy) image pair against every ROI in ``roi_zip_path``.

    Returns the list of output paths plus the bbox-adjusted ROIs (so the
    caller can build a combined ``roi.zip`` for the whole sample).
    """
    gfp_stack = _load_stack(gfp_path)
    cy_stack = _load_stack(cy_path)

    if gfp_stack.shape[1:] != cy_stack.shape[1:]:
        raise ValueError(
            f"Image pair {gfp_path.name} / {cy_path.name} has mismatched "
            f"XY dims: gfp={gfp_stack.shape}, cy={cy_stack.shape}"
        )

    rois = roifile.roiread(str(roi_zip_path))
    if isinstance(rois, roifile.ImagejRoi):
        rois = [rois]
    if not rois:
        log.warning("ROI zip %s is empty — no crops produced", roi_zip_path)
        return PairOutputs([], [], []), []

    out_gfp_dir.mkdir(parents=True, exist_ok=True)
    out_cy_dir.mkdir(parents=True, exist_ok=True)
    out_roi_dir.mkdir(parents=True, exist_ok=True)

    image_hw = gfp_stack.shape[1:]
    out_gfp: list[Path] = []
    out_cy: list[Path] = []
    out_roi_paths: list[Path] = []
    adjusted_rois: list[roifile.ImagejRoi] = []

    idx = start_index
    for roi in rois:
        bbox = _roi_bbox(roi, image_hw)
        x0, y0, x1, y1 = bbox
        mask2d, _local_poly = _local_polygon_mask(roi, bbox)

        # Crop both stacks to the bbox, then zero pixels outside the ROI
        # within the bbox (matches ImageJ "Make Inverse" + "Clear stack").
        gfp_crop = gfp_stack[:, y0:y1, x0:x1].copy()
        cy_crop = cy_stack[:, y0:y1, x0:x1].copy()
        outside = ~mask2d
        if outside.any():
            gfp_crop[:, outside] = 0
            cy_crop[:, outside] = 0

        gfp_out = out_gfp_dir / f"{gfp_prefix}{idx}.tif"
        cy_out = out_cy_dir / f"{cy_prefix}{idx}.tif"
        roi_out = out_roi_dir / f"{idx}.roi"

        # Preserve dtype of the source stack (matches the decon output, uint16).
        _save_stack(gfp_out, gfp_crop, gfp_stack.dtype)
        _save_stack(cy_out, cy_crop, cy_stack.dtype)

        adjusted = _adjust_roi_to_bbox(roi, bbox)
        adjusted.name = f"cell_{idx:04d}"
        if roi_out.exists():
            roi_out.unlink()
        adjusted.tofile(str(roi_out))

        out_gfp.append(gfp_out)
        out_cy.append(cy_out)
        out_roi_paths.append(roi_out)
        adjusted_rois.append(adjusted)
        idx += 1

    return PairOutputs(out_gfp, out_cy, out_roi_paths), adjusted_rois


def discover_sample_folders(
    root: Path,
    gfp_dirname: str,
    cy_dirname: str,
    roi_dirname: str,
) -> list[Path]:
    """Find every folder under ``root`` containing all of gfp/, cy/, roi/.

    Tolerates the user's nested ``Decon/Deconvoluted/<sample>/...`` layout
    as well as flatter ones.
    """
    root = Path(root)
    samples: list[Path] = []
    seen: set[Path] = set()
    for path in root.rglob(gfp_dirname):
        if not path.is_dir():
            continue
        sample = path.parent
        if sample in seen:
            continue
        if not (sample / cy_dirname).is_dir():
            continue
        if not (sample / roi_dirname).is_dir():
            continue
        # Skip anything already inside a Cropped/ output tree.
        if "Cropped" in sample.parts:
            continue
        samples.append(sample)
        seen.add(sample)
    return sorted(samples)


def crop_sample_folder(
    sample: Path,
    *,
    gfp_dirname: str = "gfp",
    cy_dirname: str = "cy",
    roi_dirname: str = "roi",
    gfp_prefix: str = "gfp",
    cy_prefix: str = "cy",
) -> PairOutputs:
    """Run cropping for one sample (one folder containing gfp/cy/roi)."""
    gfp_dir = sample / gfp_dirname
    cy_dir = sample / cy_dirname
    roi_dir = sample / roi_dirname

    gfp_files = _sorted_files(gfp_dir, "*.tif")
    cy_files = _sorted_files(cy_dir, "*.tif")
    roi_files = _sorted_files(roi_dir, "*.zip")

    n = min(len(gfp_files), len(cy_files), len(roi_files))
    if not (len(gfp_files) == len(cy_files) == len(roi_files)):
        log.warning(
            "%s: image/ROI count mismatch (gfp=%d, cy=%d, roi=%d) — "
            "processing the first %d only",
            sample, len(gfp_files), len(cy_files), len(roi_files), n,
        )
    if n == 0:
        log.info("%s: nothing to crop (empty channel or roi folder)", sample)
        return PairOutputs([], [], [])

    out_root = sample / "Cropped"
    out_gfp_dir = out_root / "gfp"
    out_cy_dir = out_root / "cy"
    out_roi_dir = out_root / "roi"

    all_gfp: list[Path] = []
    all_cy: list[Path] = []
    all_roi_paths: list[Path] = []
    combined_rois: list[roifile.ImagejRoi] = []

    next_idx = 1
    for gfp_path, cy_path, roi_zip in zip(gfp_files[:n], cy_files[:n], roi_files[:n]):
        log.info("Cropping pair %s + %s with %s", gfp_path.name, cy_path.name, roi_zip.name)
        outs, adjusted = crop_pair_with_rois(
            gfp_path, cy_path, roi_zip,
            out_gfp_dir, out_cy_dir, out_roi_dir,
            start_index=next_idx,
            gfp_prefix=gfp_prefix,
            cy_prefix=cy_prefix,
        )
        all_gfp.extend(outs.gfp)
        all_cy.extend(outs.cy)
        all_roi_paths.extend(outs.roi)
        combined_rois.extend(adjusted)
        next_idx += len(outs.gfp)

    # Save a combined roi.zip for the whole sample so the bg-sub stage can
    # read a single zip in image-paired order (matches the ImageJ
    # Background_Subtraction_Tailored.ijm contract).
    if combined_rois:
        combined_zip = out_root / "roi.zip"
        if combined_zip.exists():
            combined_zip.unlink()
        roifile.roiwrite(str(combined_zip), combined_rois)
        log.info("Wrote combined roi zip %s (%d ROIs)", combined_zip, len(combined_rois))

    return PairOutputs(all_gfp, all_cy, all_roi_paths)
