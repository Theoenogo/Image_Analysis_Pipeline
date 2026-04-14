"""Consolidate per-sample ``Decon/`` folders — port of MATLAB E.

After :mod:`stacking`, each per-sample folder contains a ``Decon/`` child
holding that sample's stacked ``gfp/`` and ``cy/`` outputs. This step
pulls every sample's ``Decon`` up to a single ``<main_folder>/Decon/``,
renaming each one to the originating sample folder's name so the final
structure is ``<main_folder>/Decon/<sample>/{gfp,cy}/<name>.tif``.

Matches the behavior of ``E_move_decon_folders.m``.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

log = logging.getLogger(__name__)


def consolidate_decon_folders(main_folder: Path) -> list[Path]:
    """Move every ``<sample>/Decon/`` into ``<main_folder>/Decon/<sample>/``.

    Returns the list of newly-created consolidated sample folders.
    """
    main_folder = Path(main_folder)
    if not main_folder.is_dir():
        raise NotADirectoryError(main_folder)

    top_decon = main_folder / "Decon"
    top_decon.mkdir(exist_ok=True)

    moved: list[Path] = []
    for sample_dir in sorted(p for p in main_folder.iterdir() if p.is_dir()):
        if sample_dir.name == "Decon":
            continue
        sample_decon = sample_dir / "Decon"
        if not sample_decon.is_dir():
            continue

        destination = top_decon / sample_dir.name
        if destination.exists():
            log.warning(
                "Cannot consolidate %s: %s already exists", sample_decon, destination,
            )
            continue

        shutil.move(str(sample_decon), str(destination))
        moved.append(destination)
        log.info("Moved %s -> %s", sample_decon, destination)

    return moved
