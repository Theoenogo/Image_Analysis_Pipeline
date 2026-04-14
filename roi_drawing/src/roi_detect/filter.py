"""Per-cell intensity measurement and dual-positive filtering."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class CellMeasurement:
    label: int
    area: int
    gfp_mean: float
    cy5_mean: float
    gfp_top: float  # top-percentile intensity (for punctate signal)
    cy5_top: float
    gfp_bright: float  # mean of brightest N% of pixels
    cy5_bright: float
    gfp_contrast: float  # within-cell contrast ratio (bright / median)
    cy5_contrast: float
    gfp_positive: bool
    cy5_positive: bool
    kept: bool


def _parse_metric(metric: str) -> tuple[str, float]:
    """Parse an intensity metric spec.

    Returns ``(kind, param)`` where *kind* is ``"mean"``, ``"top"``,
    ``"bright"``, or ``"contrast"`` and *param* is the numeric parameter.
    """
    metric = metric.strip().lower()
    if metric == "mean":
        return ("mean", 0.0)
    if metric == "contrast":
        return ("contrast", 10.0)  # default: top 10% vs median
    if metric.startswith("top"):
        try:
            pct = float(metric[3:])
        except ValueError:
            raise ValueError(
                f"Bad intensity-metric spec: {metric!r}. "
                "Use 'mean', 'topN', 'brightN', or 'contrast'."
            )
        if not (0 < pct <= 100):
            raise ValueError(f"Percentile must be in (0, 100], got {pct}")
        return ("top", pct)
    if metric.startswith("bright"):
        try:
            pct = float(metric[6:])
        except ValueError:
            raise ValueError(
                f"Bad intensity-metric spec: {metric!r}. "
                "Use 'mean', 'topN', 'brightN', or 'contrast'."
            )
        if not (0 < pct <= 100):
            raise ValueError(f"Percent must be in (0, 100], got {pct}")
        return ("bright", pct)
    raise ValueError(
        f"Unrecognized intensity-metric: {metric!r}. "
        "Use 'mean', 'topN' (e.g. 'top95'), 'brightN' (e.g. 'bright10'), "
        "or 'contrast'."
    )


def _bright_mean(pixels: np.ndarray, bright_pct: float) -> float:
    """Mean of the brightest *bright_pct*% of pixel values."""
    n = max(1, int(np.ceil(len(pixels) * bright_pct / 100.0)))
    # Partial sort: only need the n brightest — faster than full sort.
    top_n = np.partition(pixels, -n)[-n:]
    return float(np.mean(top_n))


def _contrast_ratio(pixels: np.ndarray, bright_pct: float) -> float:
    """Ratio of (mean of brightest N%) to (median) within a cell.

    Computed on raw pixel values (NOT background-subtracted) so the
    denominator is always a meaningful positive number. A cell with
    bright puncta on dim cytoplasm gives a ratio of ~2–4. A cell with
    uniform autofluorescence gives ~1.0–1.2.
    """
    median = float(np.median(pixels))
    if median <= 0:
        return 0.0
    bright = _bright_mean(pixels, bright_pct)
    return bright / median


def measure_cells(
    masks: np.ndarray,
    gfp: np.ndarray,
    cy5: np.ndarray,
    top_percentile: float = 95.0,
    bright_pct: float = 10.0,
) -> tuple[list[dict], float, float]:
    """Measure per-cell intensities in both channels.

    For each cell, computes four metrics:
    - **mean**, **top**, **bright**: background-subtracted intensity metrics.
    - **contrast**: ratio of (mean of brightest N%) to (median) within
      the cell, computed on raw pixel values. This is invariant to
      background level and directly measures whether the cell has
      concentrated bright signal (puncta) above its own baseline.

    Returns ``(rows, gfp_background, cy5_background)``.
    """
    from skimage.measure import regionprops

    background_mask = masks == 0
    gfp_bg = float(np.median(gfp[background_mask])) if background_mask.any() else 0.0
    cy5_bg = float(np.median(cy5[background_mask])) if background_mask.any() else 0.0

    if masks.max() == 0:
        return [], gfp_bg, cy5_bg

    props = regionprops(masks)

    rows: list[dict] = []
    for prop in props:
        coords = prop.coords  # (N, 2) array of (row, col)
        gfp_pixels = gfp[coords[:, 0], coords[:, 1]]
        cy5_pixels = cy5[coords[:, 0], coords[:, 1]]

        rows.append(
            {
                "label": int(prop.label),
                "area": int(prop.area),
                "gfp_mean": float(np.mean(gfp_pixels)) - gfp_bg,
                "cy5_mean": float(np.mean(cy5_pixels)) - cy5_bg,
                "gfp_top": float(np.percentile(gfp_pixels, top_percentile)) - gfp_bg,
                "cy5_top": float(np.percentile(cy5_pixels, top_percentile)) - cy5_bg,
                "gfp_bright": _bright_mean(gfp_pixels, bright_pct) - gfp_bg,
                "cy5_bright": _bright_mean(cy5_pixels, bright_pct) - cy5_bg,
                "gfp_contrast": _contrast_ratio(gfp_pixels, bright_pct),
                "cy5_contrast": _contrast_ratio(cy5_pixels, bright_pct),
            }
        )
    return rows, gfp_bg, cy5_bg


def resolve_threshold(spec: str | float, values: np.ndarray) -> float:
    """Turn a threshold spec into an absolute intensity cutoff.

    Accepted specs:
    - ``"all"`` or ``"none"``: accept every cell (threshold = -inf).
      Use this to disable one channel's filter entirely.
    - ``"auto"`` or ``"otsu"``: Otsu's method on the per-cell mean
      distribution. Returns 0.0 if fewer than 2 cells.
    - ``"pNN"`` (e.g. ``"p10"``, ``"p75"``): the NNth percentile of the
      per-cell means.
    - ``"zNN"`` (e.g. ``"z2"``, ``"z2.5"``): a **per-image robust
      z-score** cutoff. The threshold is ``median + z * 1.4826 * MAD``
      where MAD is the median absolute deviation of the per-cell means.
    - A float or a numeric string: used verbatim as the threshold.
    """
    if isinstance(spec, (int, float)):
        return float(spec)

    text = str(spec).strip().lower()

    if text in ("all", "none", "skip"):
        return float("-inf")

    if text in ("auto", "otsu"):
        if values.size < 2:
            return 0.0
        from skimage.filters import threshold_otsu

        try:
            return float(threshold_otsu(values))
        except ValueError:
            return 0.0

    if text.startswith("p") and text[1:].replace(".", "", 1).isdigit():
        pct = float(text[1:])
        return float(np.percentile(values, pct))

    if text.startswith("z"):
        rest = text[1:]
        # Allow a leading sign so "z-1" and "z+2" both parse.
        sign_stripped = rest[1:] if rest[:1] in ("+", "-") else rest
        if sign_stripped.replace(".", "", 1).isdigit():
            try:
                z = float(rest)
            except ValueError:
                pass
            else:
                if values.size == 0:
                    return 0.0
                median = float(np.median(values))
                mad = float(np.median(np.abs(values - median)))
                # 1.4826 scales MAD to be a consistent estimator of std
                # for a normal distribution.
                robust_sigma = 1.4826 * mad
                return median + z * robust_sigma

    try:
        return float(text)
    except ValueError as e:
        raise ValueError(
            f"Unrecognized threshold spec: {spec!r}. "
            "Use 'all', 'auto', 'pNN' (e.g. 'p10'), 'zNN' (e.g. 'z2'), "
            "or a numeric value."
        ) from e


def score_cells(
    rows: list[dict],
    gfp_threshold_spec: str | float,
    cy5_threshold_spec: str | float,
    intensity_metric: str = "mean",
) -> tuple[list[CellMeasurement], float, float]:
    """Apply dual-positive thresholds to the measured cells.

    ``intensity_metric`` selects which per-cell value the thresholds act
    on: ``"mean"``, ``"topN"``, ``"brightN"``, or ``"contrast"``.

    For ``"contrast"``, the values are dimensionless ratios (bright/median
    within each cell), so thresholds should be ratios too (e.g. ``1.5``
    or ``auto``).

    Returns ``(measurements, gfp_threshold, cy5_threshold)``.
    """
    if not rows:
        return [], 0.0, 0.0

    kind, _pct = _parse_metric(intensity_metric)
    _KEY_MAP = {"mean": "mean", "top": "top", "bright": "bright", "contrast": "contrast"}
    suffix = _KEY_MAP[kind]
    gfp_key = f"gfp_{suffix}"
    cy5_key = f"cy5_{suffix}"

    gfp_values = np.array([r[gfp_key] for r in rows], dtype=np.float64)
    cy5_values = np.array([r[cy5_key] for r in rows], dtype=np.float64)

    gfp_threshold = resolve_threshold(gfp_threshold_spec, gfp_values)
    cy5_threshold = resolve_threshold(cy5_threshold_spec, cy5_values)

    measurements: list[CellMeasurement] = []
    for r in rows:
        gfp_pos = r[gfp_key] >= gfp_threshold
        cy5_pos = r[cy5_key] >= cy5_threshold
        measurements.append(
            CellMeasurement(
                label=r["label"],
                area=r["area"],
                gfp_mean=r["gfp_mean"],
                cy5_mean=r["cy5_mean"],
                gfp_top=r["gfp_top"],
                cy5_top=r["cy5_top"],
                gfp_bright=r["gfp_bright"],
                cy5_bright=r["cy5_bright"],
                gfp_contrast=r["gfp_contrast"],
                cy5_contrast=r["cy5_contrast"],
                gfp_positive=bool(gfp_pos),
                cy5_positive=bool(cy5_pos),
                kept=bool(gfp_pos and cy5_pos),
            )
        )

    return measurements, gfp_threshold, cy5_threshold
