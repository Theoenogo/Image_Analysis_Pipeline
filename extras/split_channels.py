"""Split mixed-channel EVOS acquisitions into per-channel subfolders.

Each acquisition folder (e.g. ``con_1.2026-07-15-15-32-48``) currently holds
TIF files for both channels side-by-side.  This script moves them into
``cy{N}/`` (d1 = Cy5) and ``gfp{N}/`` (d3 = GFP) subfolders created directly
inside the parent condition folder, then removes the now-empty acquisition
folder.

Input layout (as produced by the EVOS microscope)::

    <condition>/
        <prefix>_1.<timestamp>/
            <name>d1.TIF   (Cy5 channel)
            <name>d3.TIF   (GFP channel)
            <prefix>_1.scanprotocol
        <prefix>_2.<timestamp>/
            ...

Output layout::

    <condition>/
        cy1/
            <name>d1.TIF
        gfp1/
            <name>d3.TIF
        cy2/
            <name>d2.TIF
        gfp2/
            <name>d3.TIF

Usage::

    python split_channels.py <root_dir> [--dry-run]

``root_dir`` is the top-level experiment folder (e.g. ``20260722_.../1hr/``
or the condition subfolder itself).  The script recurses through all
sub-directories automatically.
"""

from __future__ import annotations

import argparse
import logging
import re
import shutil
import sys
from pathlib import Path

log = logging.getLogger(__name__)

# Extracts the trailing integer from folder names like "con_3.2026-07-15-..." or "OA_12.2026-..."
_SAMPLE_NUM_RE = re.compile(r"_(\d+)\.")

# Maps the channel tag in the filename to its output subfolder prefix.
_CHANNEL_MAP = {
    "d1": "cy",   # Cy5
    "d3": "gfp",  # GFP
}

# Regex that matches the channel tag just before the file extension.
_CHANNEL_TAG_RE = re.compile(r"(d\d+)\.tiff?$", re.IGNORECASE)


def _is_tif(path: Path) -> bool:
    return path.suffix.lower() in (".tif", ".tiff")


def _find_acquisition_folders(root: Path) -> list[Path]:
    """Return every folder (anywhere under root) that directly contains TIF files."""
    hits: list[Path] = []
    for folder in sorted(root.rglob("*")):
        if not folder.is_dir():
            continue
        if any(_is_tif(f) for f in folder.iterdir() if f.is_file()):
            hits.append(folder)
    return hits


def split_channels(root: Path, dry_run: bool = False) -> int:
    """Walk root, split d1/d3 TIFs into channel folders in the condition directory.

    TIFs are moved from each acquisition folder into cy{N}/ and gfp{N}/ siblings
    of that acquisition folder (i.e. directly inside the condition folder).
    The acquisition folder is then removed.

    Returns the number of acquisition folders processed.
    """
    acq_folders = _find_acquisition_folders(root)

    if not acq_folders:
        log.warning("No acquisition folders with TIF files found under %s", root)
        return 0

    processed = 0
    for folder in acq_folders:
        m = _SAMPLE_NUM_RE.search(folder.name)
        if not m:
            log.warning("Cannot extract sample number from %r — skipping", folder.name)
            continue

        n = m.group(1)
        condition_dir = folder.parent  # e.g. Control/, Forskolin/, OA/

        tifs = sorted(f for f in folder.iterdir() if f.is_file() and _is_tif(f))
        if not tifs:
            continue

        # Bin each TIF by its channel tag.
        binned: dict[str, list[Path]] = {}
        unknown: list[Path] = []
        for tif in tifs:
            tag_match = _CHANNEL_TAG_RE.search(tif.name)
            if tag_match:
                tag = tag_match.group(1).lower()
                binned.setdefault(tag, []).append(tif)
            else:
                unknown.append(tif)

        rel = folder.relative_to(root)
        log.info("Processing: %s", rel)

        for tag, files in sorted(binned.items()):
            prefix = _CHANNEL_MAP.get(tag, tag)  # fall back to the raw tag if unexpected
            dest_dir = condition_dir / f"{prefix}{n}"
            log.info("  %s → %s/ (%d file(s))", tag, dest_dir.name, len(files))
            if not dry_run:
                dest_dir.mkdir(exist_ok=True)
                for f in files:
                    shutil.move(str(f), dest_dir / f.name)

        if unknown:
            log.warning(
                "  %d TIF(s) with no recognised channel tag left in place: %s",
                len(unknown),
                [f.name for f in unknown],
            )

        # Remove the now-empty acquisition folder (non-TIF files like .scanprotocol are discarded).
        if not dry_run and not unknown:
            shutil.rmtree(folder)
            log.info("  Removed: %s", folder.name)

        processed += 1

    return processed


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "root_dir",
        type=Path,
        help="Top-level experiment folder to reorganise (searched recursively)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would happen without moving any files",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show debug-level log messages",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s",
        stream=sys.stdout,
    )

    if not args.root_dir.is_dir():
        parser.error(f"Not a directory: {args.root_dir}")

    if args.dry_run:
        log.info("=== DRY RUN — no files will be moved ===")

    n = split_channels(args.root_dir, dry_run=args.dry_run)
    log.info("Done. %d acquisition folder(s) processed.", n)

    if args.dry_run:
        log.info("Re-run without --dry-run to apply changes.")


if __name__ == "__main__":
    main()
