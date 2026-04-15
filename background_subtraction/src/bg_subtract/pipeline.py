"""Per-cell background subtraction — port of ``Background_Subtraction_Tailored.ijm``.

Algorithm (matches the original ImageJ macro exactly):

1. For each ``Cropped/`` folder produced by ``roi_cropping/``, load the
   sorted gfp/ and cy/ TIFF lists and the combined ``roi.zip``.
2. For each image pair ``i`` (sorted by numeric index in the filename):
   - Select ROI ``i`` from the zip.
   - Measure the mean intensity *inside the ROI*, on the first slice only
     (matches the macro: it calls Measure once on the active slice without
     cycling through Z).
   - Compute ``subtract_value = mean * multiplier``, clamped to
     ``[100, 5000]``.
   - Subtract that constant from every pixel of every slice.
3. Save to ``<sample>/Background_Subtracted/{gfp,cy}/<filename>.tif``,
   preserving the source filename.

The output is what ``manders_mcc/`` consumes.
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


def _save_stack(path: Path, stack: np.ndarray, dtype: np.dtype) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(str(path), stack.astype(dtype, copy=False),
                     photometric="minisblack", compression=None)


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


def _roi_mask_for_image(roi: roifile.ImagejRoi, image_hw: tuple[int, int]) -> np.ndarray:
    """Rasterize the ROI into a full-image (h, w) boolean mask."""
    h, w = image_hw
    poly = _roi_polygon_xy(roi)
    rr, cc = _polygon_rasterize(poly[:, 1], poly[:, 0], shape=(h, w))
    mask = np.zeros((h, w), dtype=bool)
    mask[rr, cc] = True
    return mask


def _compute_subtract_value(
    stack: np.ndarray,
    roi: roifile.ImagejRoi,
    multiplier: float,
    *,
    floor: float = 100.0,
    ceiling: float = 5000.0,
) -> float:
    """Match ImageJ macro: mean of slice 0 inside ROI, * multiplier, clamped."""
    mask = _roi_mask_for_image(roi, stack.shape[1:])
    if not mask.any():
        return floor
    mean_intensity = float(stack[0][mask].mean())
    value = mean_intensity * multiplier
    return float(np.clip(value, floor, ceiling))


def _subtract_constant(stack: np.ndarray, value: float) -> np.ndarray:
    """Subtract ``value`` from every pixel; clip negatives to 0 (uint16-friendly)."""
    # Operate in float to allow negative intermediates, then clip at zero so
    # the cast back to uint16 doesn't wrap. Matches ImageJ's
    # ``run("Subtract...", "value=N")`` on a uint16 image.
    out = stack.astype(np.float32) - float(value)
    np.maximum(out, 0.0, out=out)
    return out


# ---------------------------------------------------------------------------
# Per-sample orchestration
# ---------------------------------------------------------------------------


@dataclass
class SampleOutputs:
    gfp: list[Path]
    cy: list[Path]


def subtract_sample_folder(
    cropped_root: Path,
    *,
    gfp_dirname: str = "gfp",
    cy_dirname: str = "cy",
    roi_zip_name: str = "roi.zip",
    gfp_multiplier: float = 1.25,
    cy_multiplier: float = 1.25,
    floor: float = 100.0,
    ceiling: float = 5000.0,
) -> SampleOutputs:
    """Run bg-subtract for one ``Cropped/`` folder.

    The output goes to a sibling ``Background_Subtracted/`` (matching the
    ImageJ macro's "shared parent" logic).
    """
    cropped_root = Path(cropped_root)
    gfp_dir = cropped_root / gfp_dirname
    cy_dir = cropped_root / cy_dirname
    roi_zip = cropped_root / roi_zip_name

    if not gfp_dir.is_dir():
        raise NotADirectoryError(gfp_dir)
    if not cy_dir.is_dir():
        raise NotADirectoryError(cy_dir)
    if not roi_zip.is_file():
        raise FileNotFoundError(roi_zip)

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

    rois = roifile.roiread(str(roi_zip))
    if isinstance(rois, roifile.ImagejRoi):
        rois = [rois]
    if len(rois) < len(gfp_files):
        raise ValueError(
            f"Not enough ROIs in {roi_zip}: rois={len(rois)}, "
            f"images={len(gfp_files)}"
        )

    out_root = cropped_root.parent / "Background_Subtracted"
    out_gfp_dir = out_root / "gfp"
    out_cy_dir = out_root / "cy"
    out_gfp_dir.mkdir(parents=True, exist_ok=True)
    out_cy_dir.mkdir(parents=True, exist_ok=True)

    out_gfp: list[Path] = []
    out_cy: list[Path] = []
    for i, (gfp_path, cy_path) in enumerate(zip(gfp_files, cy_files)):
        roi = rois[i]
        gfp_stack = _load_stack(gfp_path)
        cy_stack = _load_stack(cy_path)

        gfp_sub = _compute_subtract_value(
            gfp_stack, roi, gfp_multiplier, floor=floor, ceiling=ceiling,
        )
        cy_sub = _compute_subtract_value(
            cy_stack, roi, cy_multiplier, floor=floor, ceiling=ceiling,
        )
        log.info(
            "[%s] pair %d (%s / %s): gfp_sub=%.1f, cy_sub=%.1f",
            cropped_root.parent.name, i + 1, gfp_path.name, cy_path.name,
            gfp_sub, cy_sub,
        )

        gfp_out_arr = _subtract_constant(gfp_stack, gfp_sub)
        cy_out_arr = _subtract_constant(cy_stack, cy_sub)

        gfp_out_path = out_gfp_dir / gfp_path.name
        cy_out_path = out_cy_dir / cy_path.name
        _save_stack(gfp_out_path, gfp_out_arr, gfp_stack.dtype)
        _save_stack(cy_out_path, cy_out_arr, cy_stack.dtype)
        out_gfp.append(gfp_out_path)
        out_cy.append(cy_out_path)

    return SampleOutputs(out_gfp, out_cy)


def discover_cropped_folders(root: Path) -> list[Path]:
    """Find every ``Cropped/`` folder under ``root`` that has gfp/, cy/, roi.zip."""
    root = Path(root)
    found: list[Path] = []
    for path in root.rglob("Cropped"):
        if not path.is_dir():
            continue
        if not (path / "gfp").is_dir() or not (path / "cy").is_dir():
            continue
        if not (path / "roi.zip").is_file():
            continue
        # Skip output folders if we ever nest.
        if "Background_Subtracted" in path.parts:
            continue
        found.append(path)
    return sorted(found)


def run_pipeline(
    input_dir: Path,
    *,
    gfp_multiplier: float = 1.25,
    cy_multiplier: float = 1.25,
    floor: float = 100.0,
    ceiling: float = 5000.0,
) -> dict[Path, SampleOutputs]:
    """Walk ``input_dir`` for ``Cropped/`` folders and bg-subtract each."""
    input_dir = Path(input_dir)
    if not input_dir.is_dir():
        raise NotADirectoryError(input_dir)

    cropped = discover_cropped_folders(input_dir)
    if not cropped:
        log.warning("No Cropped/ folders found under %s", input_dir)
        return {}

    results: dict[Path, SampleOutputs] = {}
    for cropped_root in cropped:
        log.info("=== Background-subtracting: %s ===", cropped_root)
        results[cropped_root] = subtract_sample_folder(
            cropped_root,
            gfp_multiplier=gfp_multiplier,
            cy_multiplier=cy_multiplier,
            floor=floor,
            ceiling=ceiling,
        )
    return results
