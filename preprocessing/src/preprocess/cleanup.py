"""Input folder cleanup — port of MATLAB scripts A and B.

- ``rename_channel_folders`` mirrors ``A_renameFolders_v2.m``: renames
  microscope-output folders like ``gfp1.abc123`` to just ``gfp1``.
- ``delete_scan_protocol_files`` mirrors ``B_deleteScanProtocolFiles.m``:
  recursively removes ``*.scanprotocol`` sidecar files that otherwise
  interfere with folder iteration.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

log = logging.getLogger(__name__)

# The MATLAB regexes were 'gfp\d+\.[\dA-Za-z]+' etc. — match channel prefix
# followed by digits, a dot, then any alphanumeric tail. We also accept
# dashes and underscores in the tail so datetime-stamped folders like
# ``gfp2.2026-04-10-16-34-59`` are recognised and stripped.
_CHANNEL_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"^gfp\d+\.[\dA-Za-z_-]+$", re.IGNORECASE),
    re.compile(r"^cy\d+\.[\dA-Za-z_-]+$", re.IGNORECASE),
    re.compile(r"^rfp\d+\.[\dA-Za-z_-]+$", re.IGNORECASE),
)


def rename_channel_folders(main_folder: Path) -> int:
    """Strip the ``.<tail>`` suffix from channel folders anywhere in the tree.

    Walks the entire directory tree under ``main_folder`` and renames any
    folder whose name matches ``gfp\\d+.<tail>``, ``cy\\d+.<tail>``, or
    ``rfp\\d+.<tail>`` to just the part before the dot. This tolerates
    arbitrary nesting (e.g. ``<main>/experiment/replicate/sample/gfp1.abc123``).

    Returns the number of folders renamed.
    """
    main_folder = Path(main_folder)
    if not main_folder.is_dir():
        raise NotADirectoryError(main_folder)

    # Walk bottom-up so renames don't invalidate paths we still need to visit.
    to_rename: list[Path] = []
    for path in main_folder.rglob("*"):
        if not path.is_dir():
            continue
        if any(p.match(path.name) for p in _CHANNEL_PATTERNS):
            to_rename.append(path)

    renamed = 0
    for path in sorted(to_rename, key=lambda p: len(p.parts), reverse=True):
        new_name = path.name.split(".", 1)[0]
        target = path.parent / new_name
        if target.exists():
            log.warning(
                "Cannot rename %s -> %s: target already exists",
                path, target,
            )
            continue
        path.rename(target)
        renamed += 1
        log.info("Renamed %s -> %s", path, target)
    return renamed


def delete_scan_protocol_files(main_folder: Path) -> int:
    """Recursively delete every ``*.scanprotocol`` file beneath ``main_folder``.

    Returns the number of files deleted.
    """
    main_folder = Path(main_folder)
    if not main_folder.is_dir():
        raise NotADirectoryError(main_folder)

    deleted = 0
    for path in main_folder.rglob("*"):
        if path.is_file() and path.suffix.lower() == ".scanprotocol":
            path.unlink()
            deleted += 1
            log.info("Deleted %s", path)
    return deleted
