"""Rename ``rfp/`` channel folders (and their files) to ``cy/`` in place.

``roi_drawing/``, ``roi_cropping/``, ``background_subtraction/``, and
``manders_mcc/`` only recognize channel folders/files named ``gfp``/``cy``
(or ``cy5``) — none of that downstream math is fluorophore-specific, it's
just pixel-intensity thresholding, segmentation, and correlation applied to
whichever two channels get loaded. The only genuinely RFP-specific step in
the whole pipeline is the chromatic XY-offset correction in
``preprocessing/preprocess.py`` (``--gfp-offset 3 -1`` for GFP/RFP), which
happens upstream of everything this script touches.

So for a GFP/RFP dataset, once deconvolution is done, renaming ``rfp`` to
``cy`` lets every downstream stage treat it exactly like a GFP/Cy5 dataset
with no code changes. This script does that rename — both the folder name
and the leading ``rfp`` on every file inside it — since the pairing logic
in ``roi_drawing`` needs both to match (folder name for subfolder
detection, filename prefix for numeric sample-pairing).

Note the renamed dataset's CSV columns and QC plot titles will say
"Cy5"/"cy5_mean" even though the data is actually RFP — that's just a
label baked into roi_drawing's code, not a computation issue. Worth noting
in your own experiment records.

Input layout (anywhere under ``--input-dir``, typically
``Decon/Deconvoluted/<sample>/rfp/`` or ``Decon/<sample>/rfp/``)::

    <sample>/
        rfp/
            rfp01_decon.tif
            rfp02_decon.tif
            ...

Output (renamed in place)::

    <sample>/
        cy/
            cy01_decon.tif
            cy02_decon.tif
            ...

Usage::

    python rename_rfp_to_cy.py --input-dir /path/to/main_folder
    python rename_rfp_to_cy.py --input-dir /path/to/main_folder --dry-run
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

log = logging.getLogger(__name__)


def _find_rfp_folders(root: Path) -> list[Path]:
    """Find every folder anywhere under ``root`` literally named ``rfp``.

    Case-insensitive exact match only — this deliberately does *not* match
    per-slice scope-output folders like ``rfp1``/``rfp01`` (those are
    consumed earlier by the stacking step and never reach roi_drawing);
    only the stacked/deconvolved channel folder, which the pipeline always
    names exactly ``rfp``.
    """
    return sorted(
        p for p in root.rglob("*") if p.is_dir() and p.name.lower() == "rfp"
    )


def _renamed_stem(name: str) -> str | None:
    """Return ``name`` with a leading ``rfp`` replaced by ``cy``, or None."""
    if not name.lower().startswith("rfp"):
        return None
    return "cy" + name[3:]


def rename_rfp_to_cy(root: Path, dry_run: bool = False) -> list[tuple[Path, Path]]:
    """Rename every ``rfp/`` folder under ``root`` to ``cy/``, files included.

    Files are renamed first (while the parent folder still has its
    original path), then the folder itself is renamed last. Any rename
    whose target already exists is skipped with a warning rather than
    overwritten.

    Returns the list of ``(old_path, new_path)`` folder renames performed
    (or that would be performed, if ``dry_run``).
    """
    root = Path(root)
    if not root.is_dir():
        raise NotADirectoryError(root)

    rfp_folders = _find_rfp_folders(root)
    if not rfp_folders:
        log.info("No rfp/ folders found under %s", root)
        return []

    folder_renames: list[tuple[Path, Path]] = []

    for rfp_dir in rfp_folders:
        cy_dir = rfp_dir.parent / "cy"
        if cy_dir.exists():
            log.warning(
                "Skipping %s: target folder already exists: %s", rfp_dir, cy_dir,
            )
            continue

        # Rename files inside first, while the folder is still at its
        # original path.
        for file_path in sorted(rfp_dir.iterdir()):
            if not file_path.is_file():
                continue
            new_stem = _renamed_stem(file_path.name)
            if new_stem is None:
                log.warning(
                    "%s: filename doesn't start with 'rfp', leaving as-is: %s",
                    rfp_dir, file_path.name,
                )
                continue
            target = file_path.with_name(new_stem)
            if target.exists():
                log.warning(
                    "Skipping %s: target file already exists: %s",
                    file_path, target,
                )
                continue
            log.info("%s -> %s", file_path.name, target.name)
            if not dry_run:
                file_path.rename(target)

        log.info("%s -> %s", rfp_dir, cy_dir)
        if not dry_run:
            rfp_dir.rename(cy_dir)
        folder_renames.append((rfp_dir, cy_dir))

    return folder_renames


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Root folder to search for rfp/ channel folders (recursively).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be renamed without touching the filesystem.",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose logging.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stdout,
    )

    root = args.input_dir.resolve()
    if not root.is_dir():
        parser.error(f"--input-dir is not a directory: {root}")

    print(f"\n=== Rename rfp/ -> cy/ ===")
    print(f"Input folder: {root}")
    if args.dry_run:
        print("Mode: DRY RUN (no changes will be made)")
    print()

    renamed = rename_rfp_to_cy(root, dry_run=args.dry_run)

    verb = "Would rename" if args.dry_run else "Renamed"
    print(f"\n=== Done: {verb} {len(renamed)} rfp/ folder(s) ===")
    if renamed:
        for old, new in renamed:
            print(f"  {old} -> {new}")


if __name__ == "__main__":
    main()
