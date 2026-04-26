"""File IO: loading TIFFs, pairing files in a folder, saving outputs."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import tifffile


# Naming conventions we try to detect, in priority order.
# Each entry is (gfp_token, cy5_token, description). Tokens are matched
# case-insensitively as whole "words" within the filename stem, bounded by
# underscores or dashes or the start/end of the stem.
_NAMING_CONVENTIONS: list[tuple[str, str, str]] = [
    ("gfp", "cy5", "GFP / Cy5"),
    ("488", "647", "488 / 647 (wavelength)"),
    ("c1", "c2", "C1 / C2 (generic channel index)"),
    ("ch00", "ch01", "ch00 / ch01 (zero-indexed channel)"),
    ("ch01", "ch02", "ch01 / ch02 (one-indexed channel)"),
]

# Prefix-based conventions: channel name is at the START of the filename
# with no delimiter required (e.g. "gfp01_decon.tif" / "cy01_decon.tif").
_PREFIX_CONVENTIONS: list[tuple[str, str, str]] = [
    ("gfp", "cy5", "gfp/cy5 prefix"),
    ("gfp", "cy", "gfp/cy prefix"),
    ("488", "647", "488/647 prefix"),
]


@dataclass
class ImagePair:
    """A paired GFP + Cy5 acquisition.

    Either (gfp_path, cy5_path) are two separate files, or `multi_path`
    points at a single multi-channel TIFF containing both channels (with
    `gfp_index` and `cy5_index` giving their positions along the channel
    axis).
    """

    gfp_path: Path | None
    cy5_path: Path | None
    multi_path: Path | None = None
    gfp_index: int | None = None
    cy5_index: int | None = None

    @property
    def stem(self) -> str:
        """Stem to use for naming output files."""
        if self.multi_path is not None:
            return self.multi_path.stem
        # Strip the channel token from the GFP file stem to get a shared
        # base name. Try suffix conventions first, then prefix conventions.
        stem = self.gfp_path.stem  # type: ignore[union-attr]
        for gfp_tok, _cy5_tok, _desc in _NAMING_CONVENTIONS:
            stripped = _strip_token(stem, gfp_tok)
            if stripped is not None:
                return stripped
        for gfp_pfx, _cy5_pfx, _desc in _PREFIX_CONVENTIONS:
            stripped = _strip_prefix(stem, gfp_pfx)
            if stripped is not None:
                return stripped
        return stem

    def describe(self) -> str:
        if self.multi_path is not None:
            return (
                f"{self.multi_path.name} (multi-channel; "
                f"GFP=ch{self.gfp_index}, Cy5=ch{self.cy5_index})"
            )
        return f"{self.gfp_path.name}  <->  {self.cy5_path.name}"  # type: ignore[union-attr]


def _strip_token(stem: str, token: str) -> str | None:
    """Remove ``token`` from a filename stem if present as a delimited word.

    Returns the stem with the token (and one adjacent delimiter) removed,
    or None if the token is not found.
    """
    pattern = re.compile(
        rf"(?:^|(?<=[_\-\s]))({re.escape(token)})(?:$|(?=[_\-\s]))",
        re.IGNORECASE,
    )
    match = pattern.search(stem)
    if match is None:
        return None
    start, end = match.span(1)
    # Also absorb one adjacent delimiter so we don't leave dangling "__".
    if start > 0 and stem[start - 1] in "_-":
        start -= 1
    elif end < len(stem) and stem[end] in "_-":
        end += 1
    return stem[:start] + stem[end:]


def _strip_prefix(stem: str, prefix: str) -> str | None:
    """Remove a case-insensitive prefix from a filename stem.

    Handles filenames like ``gfp01_decon`` where the channel token is
    glued directly to the sample number with no delimiter.
    """
    if not stem.lower().startswith(prefix.lower()):
        return None
    rest = stem[len(prefix):]
    # Absorb a leading delimiter if present (e.g. "gfp_01" → "01").
    if rest and rest[0] in "_-":
        rest = rest[1:]
    return rest


def load_tiff(path: Path, channel_index: int | None = None) -> np.ndarray:
    """Load a TIFF as a 2-D float array.

    Handles:
    - 2-D single-plane images.
    - 3-D Z-stacks (max-projected along Z).
    - Multi-channel images (if ``channel_index`` is given, that channel
      is extracted first; the channel axis is assumed to be the smallest
      non-spatial axis).
    """
    img = tifffile.imread(str(path))
    img = np.asarray(img)

    if img.ndim == 2:
        return img.astype(np.float32)

    if channel_index is not None:
        img = _extract_channel(img, channel_index)
        if img.ndim == 2:
            return img.astype(np.float32)

    # Anything still 3-D at this point is treated as a Z-stack and
    # max-projected. Anything higher-dimensional is an error — callers
    # should be extracting channels first.
    if img.ndim == 3:
        return img.max(axis=0).astype(np.float32)

    raise ValueError(
        f"{path.name}: cannot interpret {img.ndim}-D array with shape {img.shape}. "
        "Extract a channel or Z-plane first."
    )


def _extract_channel(img: np.ndarray, channel_index: int) -> np.ndarray:
    """Extract a single channel from a multi-channel image.

    The channel axis is the smallest axis of the array (a heuristic that
    works for the common (C, H, W), (H, W, C), (C, Z, H, W), and
    (Z, C, H, W) layouts, since channel counts are typically 2–4 while
    spatial axes are much larger).
    """
    channel_axis = int(np.argmin(img.shape))
    return np.take(img, channel_index, axis=channel_axis)


def detect_channel_count(path: Path) -> int:
    """Return the smallest axis length of a TIFF — a proxy for channel count."""
    img = tifffile.imread(str(path))
    img = np.asarray(img)
    if img.ndim <= 2:
        return 1
    return int(np.min(img.shape))


# Subfolder name patterns for cross-folder pairing. Each entry is
# (gfp_folder_names, cy5_folder_names, description).
_SUBFOLDER_CONVENTIONS: list[tuple[list[str], list[str], str]] = [
    (["gfp", "488"], ["cy5", "cy", "647"], "gfp/ + cy/ subfolders"),
]


def _tiffs_in(folder: Path) -> list[Path]:
    """Return sorted TIFFs in a folder (non-recursive)."""
    if not folder.is_dir():
        return []
    return sorted(p for p in folder.iterdir() if p.suffix.lower() in (".tif", ".tiff"))


def _extract_number(stem: str) -> str | None:
    """Extract the numeric core from a filename stem for cross-folder pairing.

    Strips known channel prefixes (gfp, cy, cy5, 488, 647) and common
    suffixes (_decon, _deconvolved) to get a matching key. Purely numeric
    results are normalized to their integer value so ``"01"`` and ``"1"``
    produce the same key — files with or without leading zeros match.

    Examples:
        gfp01_decon       -> 1
        cy1_decon         -> 1
        gfp15_deconvolved -> 15
        cy5_sample03      -> sample03
        01                -> 1
    """
    s = stem.lower()
    # Strip known channel prefixes.
    for pfx in ("gfp", "cy5", "cy", "488", "647"):
        if s.startswith(pfx):
            s = s[len(pfx):]
            # Also strip a leading delimiter.
            if s and s[0] in "_-":
                s = s[1:]
            break
    # Strip common suffixes.
    for sfx in ("_deconvolved", "_decon"):
        if s.endswith(sfx):
            s = s[: -len(sfx)]
            break
    if not s:
        return None
    # Normalize purely numeric keys so "01" == "1".
    try:
        return str(int(s))
    except ValueError:
        return s


def _pair_across_folders(
    gfp_folder: Path, cy5_folder: Path
) -> list[ImagePair]:
    """Pair TIFFs across two channel-specific folders by numeric stem."""
    gfp_files = _tiffs_in(gfp_folder)
    cy5_files = _tiffs_in(cy5_folder)

    gfp_by_key: dict[str, Path] = {}
    for p in gfp_files:
        key = _extract_number(p.stem)
        if key is not None:
            gfp_by_key[key] = p

    cy5_by_key: dict[str, Path] = {}
    for p in cy5_files:
        key = _extract_number(p.stem)
        if key is not None:
            cy5_by_key[key] = p

    def _sort_key(k: str) -> tuple[int, str]:
        try:
            return (int(k), k)
        except ValueError:
            return (10**9, k)

    pairs: list[ImagePair] = []
    for key in sorted(gfp_by_key.keys(), key=_sort_key):
        if key in cy5_by_key:
            pairs.append(ImagePair(gfp_path=gfp_by_key[key], cy5_path=cy5_by_key[key]))
    return pairs


def _detect_subfolders(folder: Path) -> tuple[list[ImagePair], str | None]:
    """Check whether ``folder`` has gfp/ and cy/ (or similar) subfolders.

    Returns ``(pairs, description)`` or ``([], None)`` if no matching
    subfolder layout is found.
    """
    subdirs = {d.name.lower(): d for d in folder.iterdir() if d.is_dir()}

    for gfp_names, cy5_names, desc in _SUBFOLDER_CONVENTIONS:
        gfp_dir = None
        cy5_dir = None
        for name in gfp_names:
            if name in subdirs:
                gfp_dir = subdirs[name]
                break
        for name in cy5_names:
            if name in subdirs:
                cy5_dir = subdirs[name]
                break
        if gfp_dir is not None and cy5_dir is not None:
            pairs = _pair_across_folders(gfp_dir, cy5_dir)
            if pairs:
                return pairs, desc

    return [], None


def find_pairs(
    folder: Path,
    gfp_token: str | None = None,
    cy5_token: str | None = None,
) -> tuple[list[ImagePair], str | None]:
    """Find GFP + Cy5 image pairs in a folder.

    Returns ``(pairs, detected_convention)``. The function searches in
    order:

    1. **Subfolder layout**: if the folder contains ``gfp/`` and ``cy/``
       (or ``cy5/``, ``647/``, etc.) subdirectories, files are paired
       across them by matching the numeric part of the filename.
    2. **Same-folder layout**: files in the folder itself are paired by
       stripping a channel token from the filename stem.

    Multi-channel TIFFs (files where the smallest axis has length 2) are
    returned as single-file ImagePairs with their channel indices unset.
    """
    folder = Path(folder)

    # 1. Try subfolder layout first (gfp/ + cy/ subfolders).
    sub_pairs, sub_desc = _detect_subfolders(folder)
    if sub_pairs:
        return sub_pairs, sub_desc

    # 2. Same-folder layout.
    tiff_files = sorted(
        [p for p in folder.iterdir() if p.suffix.lower() in (".tif", ".tiff")]
    )

    # Separate multi-channel single-file images from single-channel files.
    multi_files: list[Path] = []
    single_files: list[Path] = []
    for p in tiff_files:
        try:
            ch = detect_channel_count(p)
        except Exception:
            single_files.append(p)
            continue
        if ch == 2:
            multi_files.append(p)
        else:
            single_files.append(p)

    pairs: list[ImagePair] = [
        ImagePair(gfp_path=None, cy5_path=None, multi_path=p) for p in multi_files
    ]

    if gfp_token is not None and cy5_token is not None:
        conventions = [(gfp_token, cy5_token, f"{gfp_token} / {cy5_token}")]
    else:
        conventions = _NAMING_CONVENTIONS

    detected: str | None = None
    # Try delimiter-bounded token conventions first (e.g. sample01_GFP.tif).
    for gfp_tok, cy5_tok, desc in conventions:
        matched = _pair_with_tokens(single_files, gfp_tok, cy5_tok)
        if matched:
            pairs.extend(matched)
            detected = desc
            break

    # If no suffix matches, try prefix-based conventions (e.g. gfp01_decon.tif).
    if detected is None:
        for gfp_pfx, cy5_pfx, desc in _PREFIX_CONVENTIONS:
            matched = _pair_with_prefixes(single_files, gfp_pfx, cy5_pfx)
            if matched:
                pairs.extend(matched)
                detected = desc
                break

    return pairs, detected


def _pair_with_tokens(
    files: Iterable[Path], gfp_token: str, cy5_token: str
) -> list[ImagePair]:
    """Pair a list of files using a specific (gfp_token, cy5_token) convention."""
    gfp_by_base: dict[str, Path] = {}
    cy5_by_base: dict[str, Path] = {}

    for p in files:
        stripped_gfp = _strip_token(p.stem, gfp_token)
        stripped_cy5 = _strip_token(p.stem, cy5_token)
        # A file whose stem contains BOTH tokens is ambiguous — skip it.
        if stripped_gfp is not None and stripped_cy5 is not None:
            continue
        if stripped_gfp is not None:
            gfp_by_base[stripped_gfp] = p
        elif stripped_cy5 is not None:
            cy5_by_base[stripped_cy5] = p

    pairs: list[ImagePair] = []
    for base, gfp_path in sorted(gfp_by_base.items()):
        cy5_path = cy5_by_base.get(base)
        if cy5_path is not None:
            pairs.append(ImagePair(gfp_path=gfp_path, cy5_path=cy5_path))
    return pairs


def _pair_with_prefixes(
    files: Iterable[Path], gfp_prefix: str, cy5_prefix: str
) -> list[ImagePair]:
    """Pair files where the channel token is a prefix (e.g. gfp01_decon / cy01_decon)."""
    gfp_by_base: dict[str, Path] = {}
    cy5_by_base: dict[str, Path] = {}

    for p in files:
        stripped_gfp = _strip_prefix(p.stem, gfp_prefix)
        stripped_cy5 = _strip_prefix(p.stem, cy5_prefix)
        if stripped_gfp is not None and stripped_cy5 is not None:
            continue
        if stripped_gfp is not None:
            gfp_by_base[stripped_gfp] = p
        elif stripped_cy5 is not None:
            cy5_by_base[stripped_cy5] = p

    pairs: list[ImagePair] = []
    for base, gfp_path in sorted(gfp_by_base.items()):
        cy5_path = cy5_by_base.get(base)
        if cy5_path is not None:
            pairs.append(ImagePair(gfp_path=gfp_path, cy5_path=cy5_path))
    return pairs


def load_pair(pair: ImagePair) -> tuple[np.ndarray, np.ndarray]:
    """Load both channels of an ImagePair as 2-D float32 arrays."""
    if pair.multi_path is not None:
        if pair.gfp_index is None or pair.cy5_index is None:
            raise ValueError(
                f"{pair.multi_path.name}: multi-channel file needs "
                "gfp_index and cy5_index to be set before loading."
            )
        gfp = load_tiff(pair.multi_path, channel_index=pair.gfp_index)
        cy5 = load_tiff(pair.multi_path, channel_index=pair.cy5_index)
    else:
        gfp = load_tiff(pair.gfp_path)  # type: ignore[arg-type]
        cy5 = load_tiff(pair.cy5_path)  # type: ignore[arg-type]

    if gfp.shape != cy5.shape:
        raise ValueError(
            f"{pair.describe()}: GFP and Cy5 shapes differ "
            f"({gfp.shape} vs {cy5.shape})."
        )
    return gfp, cy5


def save_masks(path: Path, masks: np.ndarray) -> None:
    """Save a labeled mask array as a 16-bit TIFF."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(str(path), masks.astype(np.uint16), compression="zlib")


def save_roi_zip(path: Path, masks: np.ndarray, labels: Iterable[int]) -> int:
    """Save a set of ImageJ polygon ROIs to a .zip file.

    Contours are extracted from the labeled mask with
    ``skimage.measure.find_contours`` (subpixel marching squares) and
    written as ImagejRoi polygon ROIs. Returns the number of ROIs written.
    """
    # Local import: skimage is heavy and only needed here.
    from skimage import measure
    import roifile

    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()

    rois: list[roifile.ImagejRoi] = []
    for label in labels:
        binary = masks == label
        if not binary.any():
            continue
        contours = measure.find_contours(binary.astype(np.float32), 0.5)
        if not contours:
            continue
        # Keep the longest contour (ignores small holes / speckles).
        contour = max(contours, key=len)
        # skimage gives (row, col) == (y, x); ImageJ ROIs want (x, y).
        points = np.column_stack([contour[:, 1], contour[:, 0]])
        roi = roifile.ImagejRoi.frompoints(points)
        roi.name = f"cell_{int(label):04d}"
        rois.append(roi)

    if rois:
        roifile.roiwrite(str(path), rois)
    return len(rois)


def save_measurements_csv(
    path: Path,
    rows: list[dict],
) -> None:
    """Write per-cell measurements to a CSV file."""
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def save_overlay(
    path: Path,
    gfp: np.ndarray,
    cy5: np.ndarray,
    masks: np.ndarray,
    kept_labels: set[int],
) -> None:
    """Write a 3-panel QC overlay PNG: GFP, Cy5, merged-with-ROIs."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage import measure

    path.parent.mkdir(parents=True, exist_ok=True)

    def _norm(img: np.ndarray) -> np.ndarray:
        lo, hi = np.percentile(img, (1.0, 99.5))
        if hi <= lo:
            return np.zeros_like(img, dtype=np.float32)
        return np.clip((img - lo) / (hi - lo), 0.0, 1.0)

    gfp_n = _norm(gfp)
    cy5_n = _norm(cy5)

    # Magenta + green merged composite (classic fluorescence false-colour).
    composite = np.stack([cy5_n, gfp_n, cy5_n], axis=-1)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes[0].imshow(gfp_n, cmap="Greens_r")
    axes[0].set_title("GFP")
    axes[1].imshow(cy5_n, cmap="magma")
    axes[1].set_title("Cy5")
    axes[2].imshow(composite)
    axes[2].set_title("Merged + ROIs")

    # Draw contours for every cell. Kept cells in green, rejected in red.
    for label in np.unique(masks):
        if label == 0:
            continue
        binary = masks == label
        contours = measure.find_contours(binary.astype(np.float32), 0.5)
        if not contours:
            continue
        contour = max(contours, key=len)
        colour = "#00ff00" if int(label) in kept_labels else "#ff3030"
        axes[2].plot(contour[:, 1], contour[:, 0], color=colour, linewidth=0.8)

    for ax in axes:
        ax.set_axis_off()

    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
