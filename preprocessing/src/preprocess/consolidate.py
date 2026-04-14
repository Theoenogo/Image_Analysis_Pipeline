"""Consolidate per-sample ``Decon/`` folders — port of MATLAB E.

After :mod:`stacking`, each folder containing channel subfolders has a
``Decon/`` child holding that sample's stacked ``gfp/`` and ``cy/``
outputs. This step pulls every such ``Decon`` up to a single
``<main_folder>/Decon/``, flattening arbitrary intermediate nesting.

For a layout like::

    <main>/experiment1/replicate1/sampleA/Decon/

the result is::

    <main>/Decon/experiment1__replicate1__sampleA/

For the simple 2-level layout (``<main>/sampleA/Decon/``) the flat name
is just ``sampleA`` — the original MATLAB behavior.

Matches the behavior of ``E_move_decon_folders.m`` but tolerates
arbitrary depth.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

log = logging.getLogger(__name__)


_FLAT_SEP = "__"


def _flatten_sample_name(sample_path: Path, main_folder: Path) -> str:
    """Turn a relative path into a flat, filesystem-safe sample identifier.

    ``<main>/expA/repA/sampleA`` -> ``"expA__repA__sampleA"``
    ``<main>/sampleA``           -> ``"sampleA"`` (unchanged for flat layouts)
    """
    rel = sample_path.relative_to(main_folder)
    return _FLAT_SEP.join(rel.parts)


def consolidate_decon_folders(main_folder: Path) -> list[Path]:
    """Move every intermediate ``Decon/`` under ``main_folder`` into ``<main>/Decon/<flat_sample>/``.

    Walks the tree, finds every ``Decon`` directory (except the top-level
    ``<main>/Decon`` itself), computes a flattened sample name from the
    parent path, and moves the folder accordingly.

    Returns the list of newly-created consolidated sample folders.
    """
    main_folder = Path(main_folder)
    if not main_folder.is_dir():
        raise NotADirectoryError(main_folder)

    top_decon = main_folder / "Decon"
    top_decon.mkdir(exist_ok=True)

    # Collect candidates first; deepest-first so we don't walk into folders
    # we're about to move.
    candidates: list[Path] = []
    for path in main_folder.rglob("Decon"):
        if not path.is_dir():
            continue
        if path == top_decon:
            continue
        # Skip anything already under the output tree.
        try:
            path.relative_to(top_decon)
            continue
        except ValueError:
            pass
        candidates.append(path)

    moved: list[Path] = []
    for sample_decon in sorted(candidates, key=lambda p: len(p.parts), reverse=True):
        sample_parent = sample_decon.parent
        flat_name = _flatten_sample_name(sample_parent, main_folder)

        destination = top_decon / flat_name
        if destination.exists():
            log.warning(
                "Cannot consolidate %s: %s already exists", sample_decon, destination,
            )
            continue

        shutil.move(str(sample_decon), str(destination))
        moved.append(destination)
        log.info("Moved %s -> %s", sample_decon, destination)

    return moved
