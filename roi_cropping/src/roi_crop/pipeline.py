"""Top-level orchestration: discover samples and crop each in turn."""

from __future__ import annotations

import logging
from pathlib import Path

from .crop import (
    PairOutputs,
    crop_sample_folder,
    discover_sample_folders,
)

log = logging.getLogger(__name__)


def run_pipeline(
    input_dir: Path,
    *,
    gfp_dirname: str = "gfp",
    cy_dirname: str = "cy",
    roi_dirname: str = "roi",
    gfp_prefix: str = "gfp",
    cy_prefix: str = "cy",
) -> dict[Path, PairOutputs]:
    """Crop every sample folder under ``input_dir``.

    Returns a mapping of sample path -> outputs. Skips sample folders that
    don't contain all three of (gfp/, cy/, roi/) at the same level.
    """
    input_dir = Path(input_dir)
    if not input_dir.is_dir():
        raise NotADirectoryError(input_dir)

    samples = discover_sample_folders(input_dir, gfp_dirname, cy_dirname, roi_dirname)
    if not samples:
        log.warning(
            "No sample folders found under %s with %s/, %s/, %s/ subfolders",
            input_dir, gfp_dirname, cy_dirname, roi_dirname,
        )
        return {}

    results: dict[Path, PairOutputs] = {}
    for sample in samples:
        log.info("=== Cropping sample: %s ===", sample)
        results[sample] = crop_sample_folder(
            sample,
            gfp_dirname=gfp_dirname,
            cy_dirname=cy_dirname,
            roi_dirname=roi_dirname,
            gfp_prefix=gfp_prefix,
            cy_prefix=cy_prefix,
        )
    return results
