"""Preprocessing for widefield microscopy images.

Ports the original MATLAB preprocessing pipeline to Python:
    1. Clean up input folders (rename scope-output folders, delete sidecar files)
    2. Stack per-slice TIFFs into multi-page stacks, applying a fixed XY offset
       to GFP to correct for chromatic / registration shift
    3. Consolidate per-sample ``Decon/`` folders into a top-level ``Decon/``
    4. Run Richardson-Lucy deconvolution against user-supplied PSFs

The original MATLAB scripts live in ``preprocessing/matlab_reference/`` for
reference.
"""

from .cleanup import delete_scan_protocol_files, rename_channel_folders
from .consolidate import consolidate_decon_folders
from .deconvolve import deconvolve_channel, richardson_lucy
from .stacking import stack_channel_folders

__all__ = [
    "rename_channel_folders",
    "delete_scan_protocol_files",
    "stack_channel_folders",
    "consolidate_decon_folders",
    "richardson_lucy",
    "deconvolve_channel",
]
