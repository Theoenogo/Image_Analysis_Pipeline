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
# followed by digits, a dot, then any alphanumeric tail.
_CHANNEL_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"^gfp\d+\.[\dA-Za-z]+$", re.IGNORECASE),
    re.compile(r"^cy\d+\.[\dA-Za-z]+$", re.IGNORECASE),
    re.compile(r"^rfp\d+\.[\dA-Za-z]+$", re.IGNORECASE),
)


def rename_channel_folders(main_folder: Path) -> int:
    """Strip the ``.<tail>`` suffix from channel folders one level down.

    Matches the MATLAB behavior: walks each first-level subdirectory of
    ``main_folder`` (the per-sample folders) and renames any child folder
    that matches ``gfp\\d+.<tail>``, ``cy\\d+.<tail>``, or ``rfp\\d+.<tail>``
    to just the part before the dot.

    Returns the number of folders renamed.
    """
    main_folder = Path(main_folder)
    if not main_folder.is_dir():
        raise NotADirectoryError(main_folder)

    renamed = 0
    for sample_dir in sorted(p for p in main_folder.iterdir() if p.is_dir()):
        for child in sorted(p for p in sample_dir.iterdir() if p.is_dir()):
            if any(p.match(child.name) for p in _CHANNEL_PATTERNS):
                new_name = child.name.split(".", 1)[0]
                target = sample_dir / new_name
                if target.exists():
                    log.warning(
                        "Cannot rename %s -> %s: target already exists",
                        child, target,
                    )
                    continue
                child.rename(target)
                renamed += 1
                log.info("Renamed %s -> %s", child, target)
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
