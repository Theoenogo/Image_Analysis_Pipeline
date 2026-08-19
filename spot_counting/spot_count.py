"""Count granule puncta in one channel and score overlap in the other.

Automates the manual scoring protocol: pick up to N puncta in the granule-marker
channel (proinsulin or mature insulin), look at the paired channel (e.g. CGA),
call overlap yes/no for each, and report a percent overlap per cell.

The perinuclear Golgi signal is **not** punctate and must never be counted. It is
excluded automatically: the reference channel is IsoData-thresholded inside the
cell, connected components at or above ``--golgi-min-area`` are called Golgi
(granule puncta are a handful of pixels; the Golgi ribbon is hundreds to
thousands), and that region is dilated so its bright rim cannot spawn false
puncta.

The two markers need different spatial rules, so ``--marker`` is **required** —
it is never inferred from the image, because both markers show perinuclear
signal plus surrounding puncta and the difference is what you intend to count:

    proinsulin  Golgi blob excluded; juxta-Golgi puncta kept (immature granules
                bud at the Golgi and move outward).
    insulin     Golgi blob excluded *and* puncta within --min-golgi-distance of
                it dropped, leaving the granules out in the processes/arms and
                along the plasma membrane.

Analysis is 2D on a max-intensity projection. Overlap is scored by nearest-
neighbour distance between independently detected spots, with --match-radius
absorbing XY drift between channels.

Input layout (any stage with ``gfp/``, ``cy/`` and ``roi.zip`` — normally
``Background_Subtracted/``, since removing diffuse haze makes puncta stand out)::

    <sample>/
        gfp/gfp01.tif  gfp02.tif  …    (one multi-page stack per cell)
        cy/cy01.tif    cy02.tif   …    (matched 1-to-1 with gfp/)
        roi.zip                         (ImageJ ROI zip, one ROI per cell)

Output (written to ``<sample>/Results/``)::

    spot_counts.csv            one row per cell
    spot_details.csv           one row per counted spot
    qc/<cell>_spots.png        overlay: Golgi zone, counted spots, other-channel spots
    qc/<cell>_spots.zip        ImageJ ROIs of the counted spots, for checking in Fiji

plus ``spot_summary.csv`` at the root of ``--input-dir`` (mean/SEM/n per sample).

Usage::

    python spot_count.py --input-dir /path/to/proinsulin_data --marker proinsulin
    python spot_count.py --input-dir /path/to/insulin_data    --marker insulin
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
from scipy import ndimage
from scipy.spatial import cKDTree
from skimage.draw import polygon as _polygon_rasterize
from skimage.feature import blob_log
from skimage.morphology import dilation, disk

log = logging.getLogger(__name__)

_NUM_RE = re.compile(r"(\d+)")

# Number of intermediate scales blob_log searches between min and max sigma.
_NUM_SIGMA = 5

# Vertices used to draw each spot's circular ROI for Fiji.
_ROI_CIRCLE_VERTICES = 16

COUNTS_COLUMNS = [
    "sample", "cell", "marker", "reference_channel",
    "spots_detected", "spots_counted", "spots_overlapping", "percent_overlap",
    "other_spots_detected", "expected_overlap_chance_pct", "chance_std_pct",
    "overlap_above_chance_pct",
    "golgi_area_px", "ref_threshold", "other_threshold",
    "min_golgi_distance_px", "max_golgi_distance_px", "match_radius_px", "seed",
]

DETAILS_COLUMNS = [
    "sample", "cell", "marker", "spot_id", "x", "y", "sigma",
    "ref_intensity", "nearest_other_distance_px", "other_intensity", "overlap",
]

SUMMARY_COLUMNS = [
    "sample", "marker", "reference_channel", "n_cells",
    "total_spots_counted", "mean_percent_overlap", "sem_percent_overlap",
    "mean_overlap_above_chance_pct", "sem_overlap_above_chance_pct",
]


@dataclass
class SpotParams:
    """Detection, exclusion and matching parameters for one run."""

    marker: str
    reference_channel: str = "gfp"
    spots_per_cell: int = 20
    seed: int = 0
    min_sigma: float = 1.0
    max_sigma: float = 3.0
    spot_threshold: float = 0.05
    golgi_min_area: int = 300
    golgi_dilate: int = 3
    match_radius: float = 3.0
    min_golgi_distance: float = 0.0
    max_golgi_distance: float = 0.0
    chance_permutations: int = 200


# ---------------------------------------------------------------------------
# File / ROI helpers
# ---------------------------------------------------------------------------

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
# Thresholding (ImageJ's Default / IsoData method)
# ---------------------------------------------------------------------------

def _isodata_from_histogram(hist: np.ndarray) -> int:
    """Core IsoData algorithm, matching ImageJ's AutoThresholder.IJDefault()."""
    n = len(hist)
    hist = np.asarray(hist, dtype=np.float64)

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
        sum2 = lo.sum()
        sum1 = (indices[min_bin:moving + 1] * lo).sum()

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


def isodata_threshold(values: np.ndarray) -> float:
    """IsoData auto-threshold for a 1-D or 2-D array of intensities."""
    values = np.asarray(values)
    if values.size == 0:
        return 0.0

    if values.dtype == np.uint8:
        hist = np.bincount(values.ravel(), minlength=256)
    elif values.dtype == np.uint16:
        hist = np.bincount(values.ravel(), minlength=65536)
    else:
        v_min, v_max = float(values.min()), float(values.max())
        if v_min == v_max:
            return v_min
        hist, _ = np.histogram(values.ravel(), bins=256, range=(v_min, v_max))
        idx = _isodata_from_histogram(hist)
        return v_min + idx * (v_max - v_min) / 256.0

    thr = _isodata_from_histogram(hist)
    if thr <= 0:
        return float(values.mean() + 2.0 * values.std())
    return float(thr)


# ---------------------------------------------------------------------------
# Image preparation
# ---------------------------------------------------------------------------

def project(stack: np.ndarray) -> np.ndarray:
    """Max-intensity-project a (Z, H, W) stack to (H, W)."""
    return stack.max(axis=0)


def normalize_in_mask(image: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Percentile-stretch to [0, 1] using percentiles measured inside the mask.

    Taking the percentiles from cell pixels only (rather than the whole frame,
    which is mostly background) is what lets one --spot-threshold work across
    images instead of needing to be retuned per image.
    """
    image = image.astype(np.float32)
    inside = image[mask] if mask.any() else image
    lo, hi = np.percentile(inside, (1.0, 99.5))
    if hi <= lo:
        return np.zeros_like(image, dtype=np.float32)
    return np.clip((image - lo) / (hi - lo), 0.0, 1.0).astype(np.float32)


def find_golgi_mask(
    projection: np.ndarray,
    cell_mask: np.ndarray,
    min_area: int,
    dilate: int,
) -> tuple[np.ndarray, float]:
    """Locate the non-punctate perinuclear region by connected-component area.

    Returns ``(golgi_mask, threshold)``. Components smaller than *min_area* are
    granule puncta and are left alone; the result is dilated by *dilate* px so
    the bright rim around the Golgi cannot produce false spots.
    """
    if not cell_mask.any():
        return np.zeros_like(cell_mask, dtype=bool), 0.0

    threshold = isodata_threshold(projection[cell_mask])
    bright = (projection >= threshold) & cell_mask

    labels, n_labels = ndimage.label(bright)
    if n_labels == 0:
        return np.zeros_like(cell_mask, dtype=bool), threshold

    areas = np.bincount(labels.ravel())
    areas[0] = 0  # background
    big_labels = np.flatnonzero(areas >= min_area)
    if big_labels.size == 0:
        return np.zeros_like(cell_mask, dtype=bool), threshold

    golgi = np.isin(labels, big_labels)
    if dilate > 0:
        golgi = dilation(golgi, disk(dilate))
    return golgi, threshold


# ---------------------------------------------------------------------------
# Spot detection & matching
# ---------------------------------------------------------------------------

def detect_spots(
    normalized: np.ndarray,
    valid_mask: np.ndarray,
    params: SpotParams,
) -> np.ndarray:
    """Detect puncta with Laplacian-of-Gaussian; keep centres inside valid_mask.

    Detection runs on the unmodified image and centres are filtered afterwards,
    rather than blanking the excluded region first — zeroing it would create a
    hard edge that the LoG filter would answer with spurious responses.

    Returns an (N, 3) array of ``(y, x, sigma)``.
    """
    blobs = blob_log(
        normalized,
        min_sigma=params.min_sigma,
        max_sigma=params.max_sigma,
        num_sigma=_NUM_SIGMA,
        threshold=params.spot_threshold,
    )
    if blobs.size == 0:
        return np.empty((0, 3), dtype=np.float64)

    yx = np.round(blobs[:, :2]).astype(int)
    h, w = normalized.shape
    in_frame = (
        (yx[:, 0] >= 0) & (yx[:, 0] < h) & (yx[:, 1] >= 0) & (yx[:, 1] < w)
    )
    blobs, yx = blobs[in_frame], yx[in_frame]
    if blobs.size == 0:
        return np.empty((0, 3), dtype=np.float64)

    keep = valid_mask[yx[:, 0], yx[:, 1]]
    return blobs[keep]


def apply_min_golgi_distance(
    spots: np.ndarray, golgi_mask: np.ndarray, min_distance: float
) -> np.ndarray:
    """Drop spots closer than *min_distance* px to the Golgi region.

    This is the mature-insulin rule: keep only granules out in the processes
    and along the plasma membrane. A no-op when no Golgi was found.
    """
    if min_distance <= 0 or spots.shape[0] == 0 or not golgi_mask.any():
        return spots

    distance = ndimage.distance_transform_edt(~golgi_mask)
    yx = np.round(spots[:, :2]).astype(int)
    return spots[distance[yx[:, 0], yx[:, 1]] >= min_distance]


def apply_max_golgi_distance(
    spots: np.ndarray, golgi_mask: np.ndarray, max_distance: float
) -> np.ndarray:
    """Keep only spots within *max_distance* px of the Golgi region.

    The opposite of apply_min_golgi_distance: for restricting a marker to
    clearly near-Golgi/immature puncta and dropping ones that have moved
    too far out into the processes to be confidently called immature. Not
    tied to a marker preset — tune by eye against the QC overlay, which
    draws this boundary as a dashed orange ring when enabled. A no-op when
    no Golgi was found.
    """
    if max_distance <= 0 or spots.shape[0] == 0 or not golgi_mask.any():
        return spots

    distance = ndimage.distance_transform_edt(~golgi_mask)
    yx = np.round(spots[:, :2]).astype(int)
    return spots[distance[yx[:, 0], yx[:, 1]] <= max_distance]


def select_spots(
    spots: np.ndarray, n_wanted: int, rng: np.random.Generator
) -> np.ndarray:
    """Randomly sample *n_wanted* spots, or return all if there are fewer.

    Random rather than brightest-first: ranking by intensity would bias toward
    large or bright granules and inflate the overlap estimate.
    """
    if spots.shape[0] <= n_wanted:
        return spots
    chosen = rng.choice(spots.shape[0], size=n_wanted, replace=False)
    return spots[np.sort(chosen)]


def match_spots(
    ref_spots: np.ndarray, other_spots: np.ndarray, radius: float
) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-neighbour distance from each reference spot to the other channel.

    Returns ``(distances, overlap)``; distance is ``inf`` where the other
    channel has no spots at all.
    """
    n = ref_spots.shape[0]
    if n == 0:
        return np.empty(0), np.empty(0, dtype=bool)
    if other_spots.shape[0] == 0:
        return np.full(n, np.inf), np.zeros(n, dtype=bool)

    tree = cKDTree(other_spots[:, :2])
    distances, _ = tree.query(ref_spots[:, :2], k=1)
    return distances, distances <= radius


def chance_level_overlap(
    ref_spots: np.ndarray,
    n_other: int,
    valid_mask: np.ndarray,
    radius: float,
    rng: np.random.Generator,
    n_permutations: int = 200,
) -> tuple[float, float]:
    """Estimate the overlap rate expected from spatial coincidence alone.

    Repeatedly scatters ``n_other`` points uniformly at random within
    ``valid_mask`` (same count and same region the real "other channel"
    spots were drawn from, so cell shape and the Golgi exclusion are
    respected) and recomputes the overlap percentage against the fixed
    ``ref_spots`` each time. This matters because a densely-detected
    "other" channel can produce a high raw percent_overlap purely from
    crowding — e.g. one real dataset here had CGA spots averaging ~9 px
    apart, which alone gives ~29% expected overlap at a 3 px match radius
    with zero true colocalization. Comparing raw percent_overlap across
    cells/datasets with different "other"-channel spot density is
    therefore misleading without this correction.

    Returns ``(mean_pct, std_pct)`` over ``n_permutations`` draws. If
    there are no reference spots or no valid area, returns ``(0.0, 0.0)``.
    """
    if ref_spots.shape[0] == 0 or n_other == 0 or not valid_mask.any():
        return 0.0, 0.0

    valid_yx = np.argwhere(valid_mask).astype(np.float64)
    n_pool = valid_yx.shape[0]
    replace = n_other > n_pool

    pcts = np.empty(n_permutations, dtype=np.float64)
    for i in range(n_permutations):
        idx = rng.choice(n_pool, size=n_other, replace=replace)
        rand_other = valid_yx[idx]
        _, overlap = match_spots(ref_spots, rand_other, radius)
        pcts[i] = 100.0 * overlap.sum() / overlap.size

    return float(pcts.mean()), float(pcts.std())


def _local_max(image: np.ndarray, y: int, x: int, radius: int) -> float:
    """Brightest pixel within *radius* of (y, x), clipped to the image."""
    h, w = image.shape
    y0, y1 = max(0, y - radius), min(h, y + radius + 1)
    x0, x1 = max(0, x - radius), min(w, x + radius + 1)
    return float(image[y0:y1, x0:x1].max())


# ---------------------------------------------------------------------------
# Per-cell analysis
# ---------------------------------------------------------------------------

def analyze_cell(
    ref_stack: np.ndarray,
    other_stack: np.ndarray,
    cell_roi: roifile.ImagejRoi,
    cell_name: str,
    sample_name: str,
    params: SpotParams,
    rng: np.random.Generator,
) -> tuple[dict, list[dict], dict]:
    """Count spots in the reference channel and score overlap in the other.

    Returns ``(counts_row, detail_rows, qc)`` where *qc* carries the arrays the
    overlay and ROI writers need.
    """
    ref_proj = project(ref_stack)
    other_proj = project(other_stack)
    hw = ref_proj.shape
    cell_mask = _roi_mask(cell_roi, hw)

    ref_norm = normalize_in_mask(ref_proj, cell_mask)
    other_norm = normalize_in_mask(other_proj, cell_mask)

    ref_golgi, ref_threshold = find_golgi_mask(
        ref_proj, cell_mask, params.golgi_min_area, params.golgi_dilate
    )
    # CGA has TGN signal too, so its Golgi is excluded on the same terms.
    other_golgi, other_threshold = find_golgi_mask(
        other_proj, cell_mask, params.golgi_min_area, params.golgi_dilate
    )

    ref_valid = cell_mask & ~ref_golgi
    other_valid = cell_mask & ~other_golgi

    ref_spots = detect_spots(ref_norm, ref_valid, params)
    ref_spots = apply_min_golgi_distance(
        ref_spots, ref_golgi, params.min_golgi_distance
    )
    ref_spots = apply_max_golgi_distance(
        ref_spots, ref_golgi, params.max_golgi_distance
    )
    other_spots = detect_spots(other_norm, other_valid, params)

    spots_detected = int(ref_spots.shape[0])
    counted = select_spots(ref_spots, params.spots_per_cell, rng)
    spots_counted = int(counted.shape[0])

    distances, overlap = match_spots(counted, other_spots, params.match_radius)
    spots_overlapping = int(overlap.sum())
    percent = (
        100.0 * spots_overlapping / spots_counted if spots_counted else 0.0
    )

    chance_mean, chance_std = chance_level_overlap(
        counted, other_spots.shape[0], other_valid, params.match_radius,
        rng, params.chance_permutations,
    )
    overlap_above_chance = percent - chance_mean

    window = max(1, int(round(params.match_radius)))
    detail_rows: list[dict] = []
    for i in range(spots_counted):
        y, x, sigma = counted[i]
        yi, xi = int(round(y)), int(round(x))
        detail_rows.append({
            "sample": sample_name,
            "cell": cell_name,
            "marker": params.marker,
            "spot_id": i + 1,
            "x": f"{x:.2f}",
            "y": f"{y:.2f}",
            "sigma": f"{sigma:.3f}",
            "ref_intensity": f"{float(ref_proj[yi, xi]):.2f}",
            "nearest_other_distance_px": (
                "" if np.isinf(distances[i]) else f"{distances[i]:.3f}"
            ),
            "other_intensity": f"{_local_max(other_proj, yi, xi, window):.2f}",
            "overlap": "yes" if overlap[i] else "no",
        })

    counts_row = {
        "sample": sample_name,
        "cell": cell_name,
        "marker": params.marker,
        "reference_channel": params.reference_channel,
        "spots_detected": spots_detected,
        "spots_counted": spots_counted,
        "spots_overlapping": spots_overlapping,
        "percent_overlap": f"{percent:.4f}",
        "other_spots_detected": int(other_spots.shape[0]),
        "expected_overlap_chance_pct": f"{chance_mean:.4f}",
        "chance_std_pct": f"{chance_std:.4f}",
        "overlap_above_chance_pct": f"{overlap_above_chance:.4f}",
        "golgi_area_px": int(ref_golgi.sum()),
        "ref_threshold": f"{ref_threshold:.2f}",
        "other_threshold": f"{other_threshold:.2f}",
        "min_golgi_distance_px": f"{params.min_golgi_distance:g}",
        "max_golgi_distance_px": f"{params.max_golgi_distance:g}",
        "match_radius_px": f"{params.match_radius:g}",
        "seed": params.seed,
    }

    qc = {
        "ref_proj": ref_proj,
        "ref_golgi": ref_golgi,
        "cell_mask": cell_mask,
        "counted": counted,
        "overlap": overlap,
        "other_spots": other_spots,
        "chance_mean": chance_mean,
        "overlap_above_chance": overlap_above_chance,
        "max_golgi_distance": params.max_golgi_distance,
    }
    return counts_row, detail_rows, qc


# ---------------------------------------------------------------------------
# QC outputs
# ---------------------------------------------------------------------------

def _spots_to_rois(
    spots: np.ndarray, overlap: np.ndarray, cell_name: str
) -> list[roifile.ImagejRoi]:
    """Draw each counted spot as a small circular polygon ROI for Fiji.

    ROI names carry the overlap call, so the ROI Manager alone tells you which
    spots were scored positive when checking against a manual count.
    """
    angles = np.linspace(0, 2 * np.pi, _ROI_CIRCLE_VERTICES, endpoint=False)
    rois: list[roifile.ImagejRoi] = []
    for i in range(spots.shape[0]):
        y, x, sigma = spots[i]
        radius = max(2.0, float(sigma) * np.sqrt(2))
        xy = np.column_stack([
            x + radius * np.cos(angles),
            y + radius * np.sin(angles),
        ]).astype(np.float32)
        tag = "ov" if overlap[i] else "no"
        rois.append(
            roifile.ImagejRoi.frompoints(xy, name=f"{cell_name}_{i + 1:02d}_{tag}")
        )
    return rois


def save_overlay(path: Path, qc: dict, cell_name: str, params: SpotParams) -> None:
    """Write the per-cell QC figure: Golgi zone, counted spots, other-channel spots."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ref_proj = qc["ref_proj"]
    counted = qc["counted"]
    overlap = qc["overlap"]

    inside = ref_proj[qc["cell_mask"]]
    vmin, vmax = (
        np.percentile(inside, (1.0, 99.5)) if inside.size else (0, 1)
    )

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(ref_proj, cmap="gray", vmin=vmin, vmax=max(vmax, vmin + 1))
    ax.contour(qc["cell_mask"], levels=[0.5], colors="white", linewidths=0.6)
    if qc["ref_golgi"].any():
        ax.contour(qc["ref_golgi"], levels=[0.5], colors="yellow", linewidths=1.0)

    max_dist = qc.get("max_golgi_distance", 0.0)
    if max_dist > 0 and qc["ref_golgi"].any():
        distance = ndimage.distance_transform_edt(~qc["ref_golgi"])
        ax.contour(distance <= max_dist, levels=[0.5], colors="orange",
                   linewidths=1.0, linestyles="dashed")

    other = qc["other_spots"]
    if other.shape[0]:
        ax.plot(other[:, 1], other[:, 0], "+", color="deepskyblue",
                markersize=5, markeredgewidth=0.8, linestyle="none")

    for i in range(counted.shape[0]):
        y, x, sigma = counted[i]
        radius = max(3.0, float(sigma) * np.sqrt(2) * 2)
        ax.add_patch(plt.Circle(
            (x, y), radius, fill=False,
            color="lime" if overlap[i] else "red", linewidth=1.2,
        ))

    n_ov = int(overlap.sum()) if overlap.size else 0
    n_tot = int(counted.shape[0])
    pct = 100.0 * n_ov / n_tot if n_tot else 0.0
    chance = qc["chance_mean"]
    above_chance = qc["overlap_above_chance"]
    legend = (
        "yellow = Golgi excluded, green = overlap, red = no overlap, "
        "blue + = other channel"
    )
    if max_dist > 0:
        legend += ", orange dashed = max Golgi distance"
    ax.set_title(
        f"{cell_name} — {params.marker} — {n_ov}/{n_tot} overlap ({pct:.1f}%, "
        f"chance ~{chance:.1f}%, {above_chance:+.1f} above chance)\n{legend}",
        fontsize=9,
    )
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Sample-level processing
# ---------------------------------------------------------------------------

def process_sample(
    sample_dir: Path,
    params: SpotParams,
    *,
    gfp_dirname: str = "gfp",
    cy_dirname: str = "cy",
    roi_zip_name: str = "roi.zip",
    roi_zip_path_override: Path | None = None,
) -> dict | None:
    """Run spot counting for one sample folder. Returns a summary dict or None.

    roi_zip_path_override
        Full path to the ROI zip when it lives outside ``sample_dir`` —
        e.g. pointing ``sample_dir`` at a ``Background_Subtracted/``
        folder while the ROI zip is still in the sibling ``Cropped/``
        folder that produced it (``background_subtraction`` never copies
        it over). Takes precedence over ``roi_zip_name``.
    """
    roi_zip = roi_zip_path_override if roi_zip_path_override is not None else sample_dir / roi_zip_name
    if not roi_zip.is_file():
        log.warning("No ROI zip at %s — skipping", roi_zip)
        return None

    gfp_dir = sample_dir / gfp_dirname
    cy_dir = sample_dir / cy_dirname
    for label, d in (("GFP", gfp_dir), ("Cy", cy_dir)):
        if not d.is_dir():
            log.warning("No %s folder in %s — skipping", label, sample_dir)
            return None

    rois = roifile.roiread(str(roi_zip))
    if isinstance(rois, roifile.ImagejRoi):
        rois = [rois]
    if not rois:
        log.warning("Empty ROI zip %s — skipping", roi_zip)
        return None

    gfp_files = _sorted_tifs(gfp_dir)
    cy_files = _sorted_tifs(cy_dir)

    if params.reference_channel == "gfp":
        ref_files, other_files = gfp_files, cy_files
    else:
        ref_files, other_files = cy_files, gfp_files

    n = min(len(ref_files), len(other_files), len(rois))
    if n == 0:
        log.warning("%s: no matched GFP/Cy/ROI triplets — skipping", sample_dir.name)
        return None
    if len(gfp_files) != len(cy_files) or len(gfp_files) != len(rois):
        log.warning(
            "%s: GFP=%d, Cy=%d, ROIs=%d — processing first %d",
            sample_dir.name, len(gfp_files), len(cy_files), len(rois), n,
        )

    out_dir = sample_dir / "Results"
    qc_dir = out_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(params.seed)
    sample_name = sample_dir.name
    counts_rows: list[dict] = []
    detail_rows: list[dict] = []

    for i in range(n):
        ref_stack = _load_stack(ref_files[i])
        other_stack = _load_stack(other_files[i])
        if ref_stack.shape[1:] != other_stack.shape[1:]:
            log.warning(
                "%s / %s: shape %s != %s — skipping",
                ref_files[i].name, other_files[i].name,
                ref_stack.shape, other_stack.shape,
            )
            continue

        cell_name = ref_files[i].stem
        counts_row, details, qc = analyze_cell(
            ref_stack, other_stack, rois[i], cell_name, sample_name, params, rng
        )
        counts_rows.append(counts_row)
        detail_rows.extend(details)

        save_overlay(qc_dir / f"{cell_name}_spots.png", qc, cell_name, params)
        spot_rois = _spots_to_rois(qc["counted"], qc["overlap"], cell_name)
        if spot_rois:
            roifile.roiwrite(
                str(qc_dir / f"{cell_name}_spots.zip"), spot_rois, mode="w"
            )

        log.info(
            "  %s: %d detected, %d counted, %d overlapping (%s%%)",
            cell_name, counts_row["spots_detected"], counts_row["spots_counted"],
            counts_row["spots_overlapping"], counts_row["percent_overlap"],
        )

    if not counts_rows:
        return None

    _write_csv(out_dir / "spot_counts.csv", COUNTS_COLUMNS, counts_rows)
    _write_csv(out_dir / "spot_details.csv", DETAILS_COLUMNS, detail_rows)

    percents = np.array([float(r["percent_overlap"]) for r in counts_rows])
    above_chance = np.array(
        [float(r["overlap_above_chance_pct"]) for r in counts_rows]
    )
    n_cells = len(percents)

    def _sem(values: np.ndarray) -> float:
        return float(values.std(ddof=1) / np.sqrt(n_cells)) if n_cells > 1 else 0.0

    sem = _sem(percents)
    above_chance_sem = _sem(above_chance)

    log.info(
        "%s: %d cells, mean overlap %.2f%% ± %.2f SEM (%.2f%% ± %.2f SEM above chance)",
        sample_name, n_cells, percents.mean(), sem,
        above_chance.mean(), above_chance_sem,
    )

    return {
        "sample": str(sample_dir),
        "marker": params.marker,
        "reference_channel": params.reference_channel,
        "n_cells": n_cells,
        "total_spots_counted": int(
            sum(r["spots_counted"] for r in counts_rows)
        ),
        "mean_percent_overlap": f"{percents.mean():.4f}",
        "sem_percent_overlap": f"{sem:.4f}",
        "mean_overlap_above_chance_pct": f"{above_chance.mean():.4f}",
        "sem_overlap_above_chance_pct": f"{above_chance_sem:.4f}",
    }


def _write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def discover_sample_folders(
    root: Path, gfp_dirname: str, cy_dirname: str, roi_zip_name: str
) -> list[Path]:
    """Find every folder under root containing gfp/, cy/ and roi.zip."""
    samples: list[Path] = []
    for path in Path(root).rglob(roi_zip_name):
        if not path.is_file():
            continue
        folder = path.parent
        if "Results" in folder.parts:
            continue
        if (folder / gfp_dirname).is_dir() and (folder / cy_dirname).is_dir():
            samples.append(folder)
    return sorted(set(samples))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

# Per-marker spatial rule. Proinsulin keeps juxta-Golgi puncta because immature
# granules bud at the Golgi; insulin drops them to leave only the mature
# granules out in the processes and at the plasma membrane.
MARKER_PRESETS = {
    "proinsulin": {"min_golgi_distance": 0.0},
    "insulin": {"min_golgi_distance": 15.0},
}


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dir", type=Path, required=True,
        help="Root folder to search for sample directories (recursively).",
    )
    parser.add_argument(
        "--marker", required=True, choices=sorted(MARKER_PRESETS),
        help="Which granule marker the reference channel holds. Selects the "
             "spatial rule and is recorded in every output row.",
    )
    parser.add_argument(
        "--reference-channel", default="gfp", choices=("gfp", "cy"),
        help="Channel holding the granule marker; spots are counted here and "
             "checked for overlap in the other channel (default: gfp).",
    )
    parser.add_argument(
        "--spots-per-cell", type=int, default=20,
        help="Maximum spots to score per cell (default: 20).",
    )
    parser.add_argument(
        "--seed", type=int, default=0,
        help="Seed for the random spot draw, for reproducibility (default: 0).",
    )
    parser.add_argument(
        "--spot-min-sigma", type=float, default=1.0,
        help="Smallest puncta scale in px (default: 1.0).",
    )
    parser.add_argument(
        "--spot-max-sigma", type=float, default=3.0,
        help="Largest puncta scale in px (default: 3.0).",
    )
    parser.add_argument(
        "--spot-threshold", type=float, default=0.05,
        help="LoG response cutoff on the normalized image. Lower to catch "
             "dimmer granules, raise to reject noise (default: 0.05).",
    )
    parser.add_argument(
        "--golgi-min-area", type=int, default=300,
        help="Connected-component area in px at or above which a bright region "
             "is Golgi rather than a granule (default: 300).",
    )
    parser.add_argument(
        "--golgi-dilate", type=int, default=3,
        help="Dilation in px around the Golgi exclusion zone (default: 3).",
    )
    parser.add_argument(
        "--match-radius", type=float, default=3.0,
        help="Overlap tolerance, absorbing XY drift between channels. In px, "
             "or nm if --pixel-size-nm is given (default: 3 px).",
    )
    parser.add_argument(
        "--min-golgi-distance", type=float, default=None,
        help="Drop spots this close to the Golgi. Defaults to the --marker "
             "preset (proinsulin: 0, insulin: 15 px). Set 0 to disable.",
    )
    parser.add_argument(
        "--max-golgi-distance", type=float, default=0.0,
        help="Keep only spots within this many px of the Golgi (0 = "
             "disabled, no cap). Opposite of --min-golgi-distance — use "
             "this to restrict a marker to clearly near-Golgi/immature "
             "puncta instead of periphery-only. Not set by --marker "
             "presets; tune by eye against the orange dashed ring in the "
             "QC overlay.",
    )
    parser.add_argument(
        "--chance-permutations", type=int, default=200,
        help="Per cell, how many times to scatter the same number of "
             "'other'-channel spots randomly within the same valid region "
             "to estimate the overlap rate expected from spatial "
             "coincidence alone (default: 200). A dense 'other' channel "
             "can produce a high percent_overlap even with zero true "
             "colocalization — this baseline is what "
             "expected_overlap_chance_pct / overlap_above_chance_pct in "
             "the output are for. Cheap to compute; raise for a tighter "
             "estimate.",
    )
    parser.add_argument(
        "--pixel-size-nm", type=float, default=None,
        help="If given, --match-radius and --min-golgi-distance are "
             "interpreted in nm and converted to px.",
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
        "--roi-zip", default="roi.zip",
        help="ROI zip filename inside each sample (default: roi.zip), or a "
             "full path to a ROI zip that lives elsewhere — e.g. the "
             "sibling Cropped/roi.zip when --input-dir points at "
             "Background_Subtracted/, since that stage never copies the "
             "zip over. A full path requires --direct (it can only apply "
             "to one sample).",
    )
    parser.add_argument(
        "--direct", action="store_true",
        help="Treat --input-dir as one sample folder instead of walking.",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose logging.",
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

    min_golgi_distance = args.min_golgi_distance
    if min_golgi_distance is None:
        min_golgi_distance = MARKER_PRESETS[args.marker]["min_golgi_distance"]

    max_golgi_distance = args.max_golgi_distance

    match_radius = args.match_radius
    if args.pixel_size_nm:
        match_radius /= args.pixel_size_nm
        min_golgi_distance /= args.pixel_size_nm
        max_golgi_distance /= args.pixel_size_nm

    params = SpotParams(
        marker=args.marker,
        reference_channel=args.reference_channel,
        spots_per_cell=args.spots_per_cell,
        seed=args.seed,
        min_sigma=args.spot_min_sigma,
        max_sigma=args.spot_max_sigma,
        spot_threshold=args.spot_threshold,
        golgi_min_area=args.golgi_min_area,
        golgi_dilate=args.golgi_dilate,
        match_radius=match_radius,
        min_golgi_distance=min_golgi_distance,
        max_golgi_distance=max_golgi_distance,
        chance_permutations=args.chance_permutations,
    )

    roi_zip_arg = Path(args.roi_zip)
    roi_zip_override: Path | None = None
    if roi_zip_arg.is_absolute() or roi_zip_arg.is_file():
        if not args.direct:
            parser.error("A full --roi-zip path requires --direct (it can only apply to one sample).")
        roi_zip_override = roi_zip_arg.resolve()
        if not roi_zip_override.is_file():
            parser.error(f"ROI zip not found: {roi_zip_override}")

    if args.direct:
        samples = [root]
    else:
        samples = discover_sample_folders(
            root, args.gfp_dirname, args.cy_dirname, args.roi_zip
        )

    print("\n=== Spot counting ===")
    print(f"Marker:            {params.marker}")
    print(f"Reference channel: {params.reference_channel}")
    print(f"Spots per cell:    {params.spots_per_cell} (seed {params.seed})")
    print(f"Match radius:      {params.match_radius:.2f} px")
    print(f"Min Golgi dist:    {params.min_golgi_distance:.2f} px")
    if params.max_golgi_distance > 0:
        print(f"Max Golgi dist:    {params.max_golgi_distance:.2f} px")
    print(f"Samples found:     {len(samples)}\n")

    if not samples:
        print("Nothing to do — no folder with gfp/, cy/ and roi.zip was found.")
        return

    summaries: list[dict] = []
    for sample in samples:
        log.info(">>> %s", sample)
        summary = process_sample(
            sample, params,
            gfp_dirname=args.gfp_dirname,
            cy_dirname=args.cy_dirname,
            roi_zip_name=args.roi_zip,
            roi_zip_path_override=roi_zip_override,
        )
        if summary:
            summaries.append(summary)

    if summaries:
        summary_path = root / "spot_summary.csv"
        _write_csv(summary_path, SUMMARY_COLUMNS, summaries)
        print(f"\nSummary written to {summary_path}")

    print(f"\n=== Done: {len(summaries)} sample(s) analyzed ===")


if __name__ == "__main__":
    main()
