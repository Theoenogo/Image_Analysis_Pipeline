"""Interactive prompts and tkinter file/folder pickers."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path


_SETTINGS_FILE = Path.home() / ".roi_detect_last.json"
_TKINTER_PROBE_CACHE: bool | None = None


@dataclass
class RunSettings:
    pan_channel: str = "gfp"  # "gfp" or "cy5"
    diameter: float | None = None
    blur_sigma: float = 2.0
    gfp_threshold: str = "p10"
    cy5_threshold: str = "auto"
    flow_threshold: float = 0.4
    cellprob_threshold: float = 0.0
    niter: int = 0  # 0 means Cellpose default (~200)
    min_size: int = 15
    segment_on: str = "pan"  # "pan" = pan-cell channel only; "both" = max of both channels
    intensity_metric: str = "mean"  # "mean" or "topN" (e.g. "top95")
    irregular_shapes: bool = False
    use_gpu: bool = True
    # Per-folder state — remembered but not core to the algorithm.
    last_folder: str = ""
    last_multi_channel_indices: dict = field(default_factory=dict)

    def save(self) -> None:
        try:
            _SETTINGS_FILE.write_text(json.dumps(asdict(self), indent=2))
        except OSError:
            pass  # Best effort; not worth failing a run over.

    @classmethod
    def load(cls) -> "RunSettings":
        if not _SETTINGS_FILE.exists():
            return cls()
        try:
            data = json.loads(_SETTINGS_FILE.read_text())
        except (OSError, json.JSONDecodeError):
            return cls()
        allowed = {f for f in cls.__dataclass_fields__}
        clean = {k: v for k, v in data.items() if k in allowed}
        return cls(**clean)


def ask(prompt: str, default: str = "") -> str:
    """Ask a free-text question with an optional default."""
    suffix = f" [{default}]" if default else ""
    answer = input(f"{prompt}{suffix}: ").strip()
    return answer or default


def ask_choice(prompt: str, choices: list[str], default: str) -> str:
    """Ask for one of a fixed set of choices (case-insensitive)."""
    lowered = [c.lower() for c in choices]
    if default.lower() not in lowered:
        raise ValueError(f"default {default!r} not in choices {choices}")
    while True:
        answer = ask(f"{prompt} ({'/'.join(choices)})", default).lower()
        if answer in lowered:
            return answer
        print(f"  Please choose one of: {', '.join(choices)}")


def ask_yes_no(prompt: str, default: bool = True) -> bool:
    default_str = "Y/n" if default else "y/N"
    while True:
        answer = input(f"{prompt} [{default_str}]: ").strip().lower()
        if not answer:
            return default
        if answer in ("y", "yes"):
            return True
        if answer in ("n", "no"):
            return False


def ask_float(prompt: str, default: float | None) -> float | None:
    default_str = "" if default is None else str(default)
    while True:
        answer = ask(prompt, default_str)
        if not answer:
            return None
        try:
            return float(answer)
        except ValueError:
            print("  Please enter a number (or press Enter for the default).")


def ask_int(prompt: str, default: int | None) -> int | None:
    default_str = "" if default is None else str(default)
    while True:
        answer = ask(prompt, default_str)
        if not answer:
            return None
        try:
            return int(answer)
        except ValueError:
            print("  Please enter an integer (or press Enter for the default).")


def _tkinter_works() -> bool:
    """Probe whether tkinter can actually create a Tk window.

    On macOS 26+ the bundled Tcl/Tk may check the host OS version against
    a minimum compiled into the library and call C abort() on mismatch.
    Python cannot catch SIGABRT via try/except, so we run the probe in a
    short-lived subprocess — if it crashes, we know to fall back to text
    prompts for the rest of the session.
    """
    global _TKINTER_PROBE_CACHE
    if _TKINTER_PROBE_CACHE is not None:
        return _TKINTER_PROBE_CACHE

    if os.environ.get("ROI_DETECT_NO_GUI"):
        _TKINTER_PROBE_CACHE = False
        return False

    try:
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                "import tkinter; r = tkinter.Tk(); r.destroy()",
            ],
            capture_output=True,
            timeout=10,
            check=False,
        )
        _TKINTER_PROBE_CACHE = result.returncode == 0
    except Exception:
        _TKINTER_PROBE_CACHE = False

    if not _TKINTER_PROBE_CACHE:
        print(
            "  (GUI file pickers unavailable on this system — using text prompts.\n"
            "   Tip: drag a file or folder from Finder into Terminal to paste its path.)"
        )

    return _TKINTER_PROBE_CACHE


def _strip_path_quotes(path_str: str) -> str:
    """Normalize a path string pasted from Finder drag-and-drop or shell."""
    s = path_str.strip()
    # Remove surrounding quotes from copy-paste.
    if len(s) >= 2 and s[0] == s[-1] and s[0] in "\"'":
        s = s[1:-1]
    # macOS Terminal escapes spaces and special chars with backslashes
    # when you drag a file in.
    for ch in " ()[]&;'\"":
        s = s.replace(f"\\{ch}", ch)
    return s


def _prompt_for_folder(title: str, initial: str = "") -> Path | None:
    """Text-mode folder prompt (used when tkinter is unavailable)."""
    print(f"\n{title}")
    while True:
        answer = ask("Folder path (drag from Finder, or blank to cancel)", initial)
        if not answer:
            return None
        path = Path(_strip_path_quotes(answer)).expanduser()
        if path.is_dir():
            return path
        print(f"  Not a folder: {path}")


def _prompt_for_file(title: str, initial: str = "") -> Path | None:
    """Text-mode file prompt (used when tkinter is unavailable)."""
    print(f"\n{title}")
    while True:
        answer = ask("File path (drag from Finder, or blank to cancel)", initial)
        if not answer:
            return None
        path = Path(_strip_path_quotes(answer)).expanduser()
        if path.is_file():
            return path
        print(f"  Not a file: {path}")


def pick_folder(title: str, initial: str = "") -> Path | None:
    """Open a tkinter folder picker, or fall back to a text prompt."""
    if not _tkinter_works():
        return _prompt_for_folder(title, initial)

    try:
        import tkinter as tk
        from tkinter import filedialog
    except ImportError:
        return _prompt_for_folder(title, initial)

    try:
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        try:
            path_str = filedialog.askdirectory(
                title=title, initialdir=initial or None
            )
        finally:
            root.destroy()
        return Path(path_str) if path_str else None
    except Exception as e:
        print(f"  ! GUI folder picker failed ({e}); using text prompt instead.")
        global _TKINTER_PROBE_CACHE
        _TKINTER_PROBE_CACHE = False
        return _prompt_for_folder(title, initial)


def pick_file(title: str, initial: str = "") -> Path | None:
    """Open a tkinter file picker, or fall back to a text prompt."""
    if not _tkinter_works():
        return _prompt_for_file(title, initial)

    try:
        import tkinter as tk
        from tkinter import filedialog
    except ImportError:
        return _prompt_for_file(title, initial)

    try:
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        try:
            path_str = filedialog.askopenfilename(
                title=title,
                initialdir=initial or None,
                filetypes=[("TIFF images", "*.tif *.tiff"), ("All files", "*.*")],
            )
        finally:
            root.destroy()
        return Path(path_str) if path_str else None
    except Exception as e:
        print(f"  ! GUI file picker failed ({e}); using text prompt instead.")
        global _TKINTER_PROBE_CACHE
        _TKINTER_PROBE_CACHE = False
        return _prompt_for_file(title, initial)


def prompt_run_settings(settings: RunSettings) -> RunSettings:
    """Prompt the user for each run parameter, using ``settings`` as defaults.

    Mutates and returns the passed-in ``RunSettings`` so the caller can
    save it to disk after a successful run.
    """
    print("\n--- Segmentation settings ---")

    settings.pan_channel = ask_choice(
        "Which channel marks ALL cells (the pan-cell channel)?",
        ["gfp", "cy5"],
        default=settings.pan_channel,
    )

    settings.irregular_shapes = ask_yes_no(
        "Are the cells non-round (e.g. INS-1 with processes, variable sizes)?",
        default=settings.irregular_shapes,
    )
    if settings.irregular_shapes:
        # Apply the preset as new defaults for the prompts below, but still
        # let the user override each one.
        from roi_detect.segment import irregular_shapes_preset
        from roi_detect.segment import SegmentationParams as _SP

        preset = irregular_shapes_preset(_SP())
        settings.blur_sigma = preset.blur_sigma
        settings.flow_threshold = preset.flow_threshold
        settings.cellprob_threshold = preset.cellprob_threshold
        settings.niter = preset.niter
        settings.segment_on = "both"
        print(
            "  Applying irregular-shapes preset: "
            f"blur_sigma={preset.blur_sigma}, "
            f"flow_threshold={preset.flow_threshold}, "
            f"cellprob_threshold={preset.cellprob_threshold}, "
            f"niter={preset.niter}, "
            f"segment_on=both."
        )
        print(
            "  Tip: for cells with processes, measure the cell BODY diameter "
            "in Fiji — not the full extent including arms."
        )

    diameter = ask_float(
        "Approximate cell BODY diameter in pixels (blank = auto-estimate)",
        settings.diameter,
    )
    settings.diameter = diameter

    settings.blur_sigma = (
        ask_float("Gaussian blur sigma for pan channel", settings.blur_sigma)
        or settings.blur_sigma
    )

    settings.segment_on = ask_choice(
        "Segment using",
        ["pan", "both"],
        default=settings.segment_on,
    )

    show_advanced = ask_yes_no(
        "Tune advanced parameters (intensity metric, Cellpose flow/cellprob/niter)?",
        default=False,
    )
    if show_advanced:
        print(
            "\n  Intensity metric: 'mean' (default), 'contrast' (best for "
            "punctate signal),\n  'brightN' (e.g. 'bright10'), or 'topN' "
            "(e.g. 'top95')."
        )
        settings.intensity_metric = ask(
            "  Intensity metric",
            settings.intensity_metric,
        )
        settings.flow_threshold = (
            ask_float(
                "flow_threshold (higher = more irregular shapes allowed, 0.0-1.0)",
                settings.flow_threshold,
            )
            or settings.flow_threshold
        )
        settings.cellprob_threshold = (
            ask_float(
                "cellprob_threshold (lower = catches dimmer cell edges, -6 to 6)",
                settings.cellprob_threshold,
            )
            if settings.cellprob_threshold is not None
            else 0.0
        )
        # ask_float returns None on blank input; coerce back to the default.
        if settings.cellprob_threshold is None:
            settings.cellprob_threshold = 0.0
        niter = ask_int(
            "niter (0=default ~200; use 2000+ for long thin processes)",
            settings.niter,
        )
        settings.niter = niter if niter is not None else 0
        min_size = ask_int(
            "min_size in pixels (drop detections smaller than this)",
            settings.min_size,
        )
        settings.min_size = min_size if min_size is not None else 15

    # Threshold defaults depend on the intensity metric and role.
    # For "contrast" mode, values are ratios (1.0 = no contrast,
    # 2.0 = bright pixels are 2x the median). For intensity modes,
    # values are background-subtracted intensities.
    is_contrast = settings.intensity_metric.strip().lower() == "contrast"

    if is_contrast:
        # Contrast metric: selective channel needs a ratio threshold,
        # pan channel can be skipped (all cells are positive).
        if settings.pan_channel == "gfp":
            gfp_thr_default = "all"
            cy5_thr_default = "1.5"
        else:
            gfp_thr_default = "1.5"
            cy5_thr_default = "all"
    elif settings.irregular_shapes:
        # Calibrated defaults for expanded masks (cellprob_threshold <= -3).
        if settings.pan_channel == "gfp":
            gfp_thr_default = "100"
            cy5_thr_default = "150"
        else:
            gfp_thr_default = "150"
            cy5_thr_default = "100"
    else:
        # Standard defaults: pan channel low-bar, selective Otsu.
        if settings.pan_channel == "gfp":
            gfp_thr_default = "p10"
            cy5_thr_default = "auto"
        else:
            gfp_thr_default = "auto"
            cy5_thr_default = "p10"

    pan_label = settings.pan_channel.upper()
    selective = "Cy5" if settings.pan_channel == "gfp" else "GFP"

    if is_contrast:
        print(
            "\nUsing 'contrast' metric — thresholds are ratios (bright/median "
            "within each cell).\n  Typical values: 1.3 (permissive) to 2.0 "
            "(strict). 'all' skips the channel."
        )
    else:
        print(
            "\nThreshold specs: 'auto' = Otsu, 'pNN' = NNth percentile "
            "(e.g. 'p10'), 'zNN' = robust z-score (e.g. 'z2'), 'all' = skip "
            "this channel (accept all cells), or a numeric value."
        )
    print(
        f"  Defaults: {selective} (selective) = {gfp_thr_default if settings.pan_channel == 'cy5' else cy5_thr_default}, "
        f"{pan_label} (pan-cell) = {gfp_thr_default if settings.pan_channel == 'gfp' else cy5_thr_default}."
    )
    settings.gfp_threshold = ask(
        "GFP threshold spec"
        + (" (selective)" if settings.pan_channel == "cy5" else " (pan-cell)"),
        gfp_thr_default,
    )
    settings.cy5_threshold = ask(
        "Cy5 threshold spec"
        + (" (selective)" if settings.pan_channel == "gfp" else " (pan-cell)"),
        cy5_thr_default,
    )

    return settings
