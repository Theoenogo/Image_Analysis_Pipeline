"""Public API for the background_subtraction package."""

from __future__ import annotations

from .pipeline import (
    discover_cropped_folders,
    run_pipeline,
    subtract_sample_folder,
)

__all__ = [
    "discover_cropped_folders",
    "run_pipeline",
    "subtract_sample_folder",
]
