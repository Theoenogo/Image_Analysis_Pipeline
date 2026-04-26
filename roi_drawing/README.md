# ROI-Drawing

Automatically detect cells that are positive in **both** the GFP and Cy5
channels of a fluorescence microscopy image pair, draw ROIs around each
cell's full membrane, and export them as ImageJ-compatible ROI `.zip`
files.

Built around [Cellpose](https://github.com/mouseland/cellpose), which
handles punctate cytoplasmic signal where classical threshold + watershed
pipelines break down.

## What it does

1. Loads a pan-cell channel (GFP or Cy5, whichever one has signal in every
   cell) and a selective channel (the other one).
2. Lightly Gaussian-blurs the pan-cell channel to bridge punctate signal
   into continuous cell bodies.
3. Runs Cellpose (`cyto3` model) to segment every cell.
4. Measures background-subtracted mean intensity of each cell in both
   channels and keeps only the ones that are bright in **both**.
5. Writes out the results for review in Fiji.

Output files per image pair:

- `<stem>_rois.zip` — ImageJ ROI set with all kept cells (open in
  Fiji's ROI Manager).
- `roi_original/<NN>.zip` — numbered copy of the ROI zip in a shared
  `roi_original/` folder: `01.zip` for image pair 1, `02.zip` for pair 2,
  etc. The next pipeline stage expects a user-edited sibling folder
  named `roi/` (produced by the ImageJ
  `Color_Merge_Automated_PreloadROIs_Adjust` macro — run manually to
  add/remove/refine ROIs), so keeping the automated output under
  `roi_original/` preserves it verbatim.
- `<stem>_masks.tif` — labeled mask of kept cells only.
- `<stem>_overlay.png` — 3-panel QC figure (GFP / Cy5 / merged with ROIs).
- `<stem>_measurements.csv` — per-cell intensities (mean, top-percentile,
  bright, and contrast), positive/negative calls, and the thresholds used.

## Install

```bash
pip install -r requirements.txt
```

First run will download the Cellpose `cyto3` model weights (~25 MB).

**Platform notes:**
- **Linux**: you may need `sudo apt install python3-tk` for the
  folder/file pickers. The script falls back to text prompts if
  tkinter isn't available.
- **Windows**: tkinter is included with the standard Python installer.
  If you used the Microsoft Store Python, tkinter should already be
  available. All file paths use `pathlib` internally, so backslash
  paths work as expected.
- **macOS**: tkinter ships with the python.org installer and
  Homebrew Python.

A GPU is optional but speeds Cellpose up considerably. Pass `--no-gpu` to
force CPU mode.

## Usage

**Interactive (the common case — no paths to type):**

```bash
python roi_detect.py
```

A folder picker pops up. The script scans the folder, auto-detects GFP/Cy5
pairs, asks you a few parameter questions, and processes everything. Your
answers are saved to `~/.roi_detect_last.json` and offered as defaults on
the next run — so once you've dialed in settings for your microscope, you
can just press Enter through the prompts.

Pick "Single" at the mode prompt to process one image pair at a time for
parameter tuning.

**Non-interactive batch:**

```bash
python roi_detect.py \
    --batch path/to/folder \
    --pan-channel gfp \
    --diameter 30
```

**Non-interactive single-pair:**

```bash
python roi_detect.py \
    --gfp path/to/sample_GFP.tif \
    --cy5 path/to/sample_Cy5.tif \
    --pan-channel gfp \
    --diameter 30
```

## File naming conventions

### Separate folders (auto-detected)

If the selected folder contains `gfp/` and `cy/` (or `cy5/`, `647/`)
subfolders, files are paired across them by matching the numeric part
of the filename:

```
Deconvoluted/          ← point the batch at this folder
  gfp/
    gfp01_decon.tif
    gfp02_decon.tif
  cy/
    cy01_decon.tif
    cy02_decon.tif
```

The script strips the channel prefix (`gfp`, `cy`, `cy5`, `488`, `647`)
and common suffixes (`_decon`, `_deconvolved`) to get a matching key
(e.g. `01`, `02`).

### Same folder

When all TIFFs are in a single folder, pairs are matched by stem. These
suffixes are auto-detected case-insensitively:

- `<stem>_GFP.tif` ↔ `<stem>_Cy5.tif`
- `<stem>_488.tif` ↔ `<stem>_647.tif`
- `<stem>_C1.tif` ↔ `<stem>_C2.tif`
- `<stem>_ch00.tif` ↔ `<stem>_ch01.tif`
- `<stem>_ch01.tif` ↔ `<stem>_ch02.tif`

Multi-channel TIFFs (a single file containing both channels) are also
supported — the script asks once which channel index is GFP vs Cy5.

## Threshold specs

`--gfp-threshold` and `--cy5-threshold` accept:

- `auto` — Otsu's method on the per-cell mean distribution. Good default
  for the **selective** channel when the positive and negative populations
  are cleanly bimodal.
- `pNN` — the NNth percentile of per-cell means, e.g. `p10`, `p75`. Good
  default for the **pan-cell** channel — use a low percentile like `p10`
  to mostly reject segmentation noise.
- `zNN` — a **robust per-image z-score** cutoff, e.g. `z2`, `z2.5`. The
  threshold is `median + z * 1.4826 * MAD` where MAD is the median
  absolute deviation of the per-cell means in *this* image. Use this
  when a fixed absolute number works on most images but fails on a few
  because background level or spread drifts image-to-image: a cell at
  intensity 107 in an image where everyone else is 20–40 has a large
  z-score and gets kept, even though it would fall below a fixed
  `130` cutoff calibrated on a brighter image. Typical values are
  `z2` (permissive) to `z3` (strict).
- `all` — accept every cell (disable this channel's filter). Use this
  to select on one channel only, e.g. `--cy5-threshold all` to skip
  the Cy5 filter and select purely on GFP.
- A number — an absolute intensity cutoff. For `contrast` metric, this
  is a ratio (e.g. `1.5`); for other metrics, it's a background-subtracted
  intensity value.

## Intensity metric

`--intensity-metric` controls which per-cell value the thresholds act on:

- `mean` (default) — the average pixel intensity within the cell mask.
  Works well when signal is distributed uniformly across the cell.
- `contrast` — **ratio of (mean of brightest 10% of pixels) to (median
  pixel)** within each cell. **Best for punctate signal with
  autofluorescence.** This metric is background-invariant: a cell with
  bright puncta on dim cytoplasm gives a ratio of ~2–4, while a cell
  with uniform autofluorescence gives ~1.0–1.2, regardless of the
  overall image brightness. Threshold with a number (e.g. `1.5`), `auto`,
  or `z2`.
- `brightN` (e.g. `bright10`) — mean of the brightest N% of pixels.
- `topN` (e.g. `top95`) — Nth percentile of pixel intensities.

Example:

```bash
python roi_detect.py --batch ./imgs --pan-channel gfp \
    --intensity-metric contrast --gfp-threshold 1.5 --cy5-threshold all
```

To skip one channel's filter entirely, set its threshold to `all`:

```bash
--cy5-threshold all    # select on GFP only
--gfp-threshold all    # select on Cy5 only
```

All metric values (mean, top, bright, contrast) are always written to
the measurements CSV regardless of which metric is used for thresholding.

## Tuning tips

- **Cell diameter** is the most important parameter. Measure a typical
  cell in Fiji (line tool + `Measure`) and pass it as `--diameter`. The
  auto-estimate is unreliable on punctate data.
- If ROIs leave gaps in the middle of cells, increase `--blur-sigma`
  (e.g. from `2.0` to `3.0` or `4.0`) to bridge more of the puncta.
- If adjacent cells are merged into one ROI, decrease `--blur-sigma`, or
  lower `--diameter`, or raise `--cellprob-threshold`.
- If too many real double-positive cells are rejected, loosen the
  selective-channel threshold (e.g. `auto` → `p50` → `p25`).
- **Punctate signal** — if `z2` is too permissive and `z2.5` too strict,
  switch to `--intensity-metric bright10` (mean of brightest 10% of
  pixels in each cell). This gives a cleaner separation than tuning the
  z-score because it measures what matters — whether the cell has
  concentrated bright spots — without being fooled by hot pixels the
  way `top95` can be.

## Non-round cells with processes (INS-1 832/13, neurons, etc.)

The defaults above are tuned for roughly round cells. If your cells have
protrusions, arms, or variable sizes — for example **INS-1 832/13**
pancreatic β-cells, primary neurons, or any cytoskeletally-active line —
pass `--irregular-shapes` to switch on a preset:

```bash
python roi_detect.py --batch path/to/folder --pan-channel gfp \
    --diameter 22 --irregular-shapes
```

The preset sets:

- `blur_sigma=1.0` — lighter blur so thin processes are preserved.
- `flow_threshold=0.6` — allows more irregular shapes through Cellpose's
  shape filter.
- `cellprob_threshold=-3.0` — more permissive pixel inclusion so dim
  process tips don't get cut off.
- `niter=4000` — runs Cellpose's flow dynamics for many more iterations
  so the flow field propagates all the way down long processes. This is
  the single most important knob for cells with arms.
- `segment_on=both` — segments on the element-wise max of both channels
  so ROIs cover the full extent of signal in either channel, not just
  the pan-cell channel.

**Important — measure the cell BODY diameter**, not the total extent
including arms. INS-1 832/13 cell bodies are typically around 15–25 μm
(convert to pixels using your pixel size). Cellpose rescales internally
based on this, so overestimating because of long processes will split
cells up.

You can override any individual value from the preset with explicit
CLI flags, e.g.:

```bash
python roi_detect.py --batch ./imgs --pan-channel gfp --diameter 22 \
    --irregular-shapes --niter 3000 --cellprob-threshold -2.0
```

The interactive prompt also asks "Are the cells non-round?" and offers
the same preset followed by per-parameter overrides.

### Tuning checklist for irregular cells

- **Processes get cut off** → increase `--niter` (e.g. `3000`, `5000`)
  and/or lower `--cellprob-threshold` (e.g. `-2.0`).
- **Cells with arms are split into multiple ROIs** → increase `--niter`
  (the flow didn't reach down the process). Also try raising
  `--flow-threshold` to `0.8`.
- **Thin processes disappear entirely** → lower `--blur-sigma` to `0.5`
  or `0`.
- **Variable cell sizes** → set `--diameter` to the *median* body
  diameter. Cellpose rescales and handles 2–3x size variation well; for
  extreme size ranges you may need to run twice with different diameters
  and merge results manually.
- **Too much debris** → raise `--min-size` (e.g. to `50` or `100`).

## Reviewing results in Fiji

1. `File > Open...` the original image (e.g. the GFP TIFF).
2. Drag `<stem>_rois.zip` onto the Fiji window — the ROIs appear in the
   ROI Manager.
3. Use the ROI Manager to show, rename, delete, or add ROIs, then
   `More >> Save...` to re-export.
