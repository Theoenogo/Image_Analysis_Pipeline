"""Stack per-slice TIFFs into multi-page stacks — port of MATLAB C2/D2.

Each per-sample folder under ``main_folder`` contains channel subfolders
(``gfp1``, ``cy1``, etc.) holding individual TIFF slices. This module
reads those slices, optionally applies a fixed XY offset (to correct for
GFP/Cy5 channel registration), writes a multi-page uint16 TIFF named after
the channel folder, and moves it to ``<sample>/Decon/<channel_group>/``.

Matches the MATLAB behavior from
``C2_green_stack_offset_Cy5_recursive.m`` (gfp, offset (+5, -2)) and
``D2_cy5_stack_no_offset_recursive.m`` (cy, no offset).
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np
import tifffile

log = logging.getLogger(__name__)


# Default XY offset applied to GFP in the MATLAB pipeline.
# Positive x shifts the image right; positive y shifts down.
DEFAULT_GFP_OFFSET: tuple[int, int] = (5, -2)

_CHANNEL_NUM_RE = re.compile(r"^([A-Za-z]+)(\d+)$")


def _apply_xy_offset(image: np.ndarray, xy_offset: tuple[int, int]) -> np.ndarray:
    """Shift a 2D image by ``(x, y)`` pixels, zero-padding the exposed edge.

    Matches MATLAB's ``imwarp`` with an ``affine2d`` of translation
    ``[xOffset yOffset 1]`` and ``'OutputView', imref2d(size(...))``:
    the output has the same shape, data shifted, and revealed borders
    filled with zeros.
    """
    dx, dy = xy_offset
    if dx == 0 and dy == 0:
        return image.copy()

    shifted = np.zeros_like(image)
    src_y_start = max(0, -dy)
    src_y_end = image.shape[0] - max(0, dy)
    src_x_start = max(0, -dx)
    src_x_end = image.shape[1] - max(0, dx)

    dst_y_start = max(0, dy)
    dst_y_end = dst_y_start + (src_y_end - src_y_start)
    dst_x_start = max(0, dx)
    dst_x_end = dst_x_start + (src_x_end - src_x_start)

    shifted[dst_y_start:dst_y_end, dst_x_start:dst_x_end] = image[
        src_y_start:src_y_end, src_x_start:src_x_end
    ]
    return shifted


def _read_slices(folder: Path) -> tuple[list[np.ndarray], list[Path]]:
    """Read every single-page ``.tif``/``.tiff`` in ``folder``, sorted by name.

    Skips the possibly-present multipage output file that has the same
    name as the folder (``<folder>/<folder_name>.tif``) so re-runs don't
    re-ingest their own output.
    """
    candidates = sorted(
        [p for p in folder.iterdir() if p.suffix.lower() in (".tif", ".tiff")],
        key=lambda p: p.name.lower(),
    )
    # Exclude the folder-named output file if it already exists (re-run safety).
    # Check both the raw folder name and the zero-padded variant.
    exclude = {f"{folder.name}.tif"}
    m = _CHANNEL_NUM_RE.match(folder.name)
    if m:
        prefix, digits = m.group(1), m.group(2)
        exclude.add(f"{prefix}{int(digits):02d}.tif")
    candidates = [p for p in candidates if p.name not in exclude]

    slices: list[np.ndarray] = []
    used: list[Path] = []
    for path in candidates:
        try:
            arr = tifffile.imread(str(path))
        except Exception as err:
            log.warning("Skipping unreadable TIFF %s: %s", path, err)
            continue
        # Skip multi-page files accidentally present in the per-slice folder.
        if arr.ndim > 2:
            log.warning(
                "Skipping %s: unexpected multi-page TIFF (ndim=%d) in per-slice folder",
                path, arr.ndim,
            )
            continue
        slices.append(arr)
        used.append(path)
    return slices, used


def _save_stack_uint16(path: Path, slices: list[np.ndarray]) -> None:
    """Save a list of 2D arrays as a multi-page uint16 TIFF."""
    stack = np.stack([np.asarray(s).astype(np.uint16) for s in slices], axis=0)
    tifffile.imwrite(str(path), stack, imagej=True, photometric="minisblack")


def _process_channel_folder(
    channel_folder: Path,
    decon_group: str,
    xy_offset: tuple[int, int] | None,
) -> Path | None:
    """Stack slices in ``channel_folder``, apply offset if given, move to Decon.

    Returns the path to the generated stack inside ``<parent>/Decon/<group>/``,
    or ``None`` if no slices were found.
    """
    slices, _used = _read_slices(channel_folder)
    if not slices:
        log.info("No TIFF slices found in %s — skipping", channel_folder)
        return None

    if xy_offset is not None:
        slices = [_apply_xy_offset(img, xy_offset) for img in slices]

    # Zero-pad the channel number in the output filename so gfp1 → gfp01.tif.
    folder_name = channel_folder.name
    m = _CHANNEL_NUM_RE.match(folder_name)
    if m:
        prefix, digits = m.group(1), m.group(2)
        if len(digits) < 2:
            folder_name = f"{prefix}{int(digits):02d}"
    output_name = f"{folder_name}.tif"
    temp_output = channel_folder / output_name
    _save_stack_uint16(temp_output, slices)

    decon_dir = channel_folder.parent / "Decon" / decon_group
    decon_dir.mkdir(parents=True, exist_ok=True)
    final_output = decon_dir / output_name
    if final_output.exists():
        final_output.unlink()
    temp_output.replace(final_output)

    log.info("Stacked %d slice(s) from %s -> %s", len(slices), channel_folder, final_output)
    return final_output


def stack_channel_folders(
    main_folder: Path,
    channel_prefix: str,
    decon_group: str,
    xy_offset: tuple[int, int] | None = None,
) -> list[Path]:
    """Stack every ``<channel_prefix>*`` folder found anywhere under ``main_folder``.

    Walks the entire directory tree and stacks every folder whose name starts
    with ``channel_prefix`` (case-insensitive). This tolerates arbitrary
    nesting (e.g. ``<main>/experiment/replicate/sample/gfp1/``). Folders
    under any existing ``<main>/Decon/`` subtree are skipped so re-runs
    don't re-ingest their own output.

    Output is written next to the sample (i.e. as a sibling of the channel
    folder under ``<sample>/Decon/<decon_group>/<folder_name>.tif``), which
    the consolidate step later moves up to ``<main>/Decon/``.

    Parameters
    ----------
    main_folder
        The top-level folder to walk.
    channel_prefix
        Case-insensitive folder-name prefix to match (e.g. ``"gfp"``, ``"cy"``).
    decon_group
        Subfolder name under ``Decon/`` to place the stacked output in
        (e.g. ``"gfp"``, ``"cy"``).
    xy_offset
        Optional ``(x, y)`` pixel shift applied per slice before stacking.
        Pass ``None`` for no shift.

    Returns the list of output stack paths.
    """
    main_folder = Path(main_folder)
    if not main_folder.is_dir():
        raise NotADirectoryError(main_folder)

    prefix = channel_prefix.lower()
    top_decon = main_folder / "Decon"
    outputs: list[Path] = []

    channel_folders: list[Path] = []
    for path in main_folder.rglob("*"):
        if not path.is_dir():
            continue
        # Skip anything under the top-level <main>/Decon/ output tree.
        try:
            path.relative_to(top_decon)
            continue
        except ValueError:
            pass
        # Skip per-sample intermediate Decon/ folders we may have already created.
        if "Decon" in path.parts[len(main_folder.parts):]:
            continue
        if path.name.lower().startswith(prefix):
            channel_folders.append(path)

    for channel_folder in sorted(channel_folders):
        out = _process_channel_folder(channel_folder, decon_group, xy_offset)
        if out is not None:
            outputs.append(out)

    return outputs
