# Extras

Standalone analysis scripts that are useful but not part of the main
pipeline (preprocessing → ROI drawing → cropping → background
subtraction → colocalization). Each script is self-contained and can
be run independently.

## Scripts

### `split_channels.py` — EVOS channel folder reorganization

Reorganizes mixed-channel EVOS acquisition folders into the per-channel
subfolder structure expected by the rest of the pipeline.

The EVOS microscope saves both channels into a single timestamped
folder (e.g. `con_1.2026-07-15-15-32-48/`) with filenames ending in
`d1` (Cy5) or `d3` (GFP). This script moves those files into
`cy{N}/` and `gfp{N}/` subfolders created directly inside the parent
condition folder, then removes the now-empty acquisition folder.

**Run this before `preprocessing/` on any EVOS dataset.**

**Input layout** (EVOS output):
```
<condition>/
    con_1.2026-07-15-15-32-48/
        <name>d1.TIF   ← Cy5
        <name>d3.TIF   ← GFP
        con_1.scanprotocol
    con_2.2026-07-15-15-45-00/
        ...
```

**Output layout** (ready for preprocessing):
```
<condition>/
    cy1/
        <name>d1.TIF
    gfp1/
        <name>d3.TIF
    cy2/
        <name>d1.TIF
    gfp2/
        <name>d3.TIF
```

**Usage:**

```bash
python split_channels.py /path/to/experiment_folder
```

Preview what would happen without moving any files:

```bash
python split_channels.py /path/to/experiment_folder --dry-run
```

**Options:**

| Flag | Purpose |
|------|---------|
| `--dry-run` | Print what would happen without moving any files |
| `-v / --verbose` | Show debug-level log messages |

The script recurses through all subdirectories automatically, so you
can point it at the top-level experiment folder and it will find every
acquisition folder underneath.

---

### `measure_roi_signal.py` — Per-slice ROI signal measurement

Measures the mean intensity, integrated density, area, min, max, and
standard deviation of pixels inside each cell's ROI on every Z slice.
Writes a single combined CSV per sample with all cells and channels.

Python port of the ImageJ macro
[`Measure_ROI_Signal_Per_Slice.ijm`](./imagej_reference/Measure_ROI_Signal_Per_Slice.ijm).
The macro wrote one CSV per cell; this script consolidates everything
into one CSV.

**Input:** any folder structure containing `{gfp,cy}/` channel
subfolders and a `roi.zip` — works with `Cropped/` output from
`roi_cropping/` or `Background_Subtracted/` output from
`background_subtraction/`.

**Output:** `<sample>/Results/signal_measurements.csv`

**CSV columns:**

| Column | Description |
|--------|-------------|
| `cell` | Cell identifier (filename stem, e.g. `gfp01`) |
| `channel` | Channel folder name (`gfp` or `cy`) |
| `image_file` | Source TIFF filename |
| `slice` | Z-slice index (1-based) |
| `area` | Number of pixels inside the ROI |
| `mean` | Mean pixel intensity inside the ROI |
| `std_dev` | Standard deviation of pixel intensities |
| `min` | Minimum pixel intensity inside the ROI |
| `max` | Maximum pixel intensity inside the ROI |
| `integrated_density` | Sum of all pixel intensities inside the ROI |

**Usage:**

**Auto-discovery mode** — point at a root folder, the script finds all
samples:

```bash
python measure_roi_signal.py --input-dir /path/to/main_folder
```

To measure only one channel:

```bash
python measure_roi_signal.py --input-dir /path/to/main_folder --channels cy
```

**Direct mode** — point at a single channel folder and give the full
path to the ROI zip (useful when the ROI zip isn't in the same
directory as the channel TIFFs):

```bash
python measure_roi_signal.py \
    --input-dir /path/to/Background_Subtracted/cy \
    --roi-zip /path/to/Cropped/roi.zip
```

**Common options:**

| Flag | Default | Purpose |
|------|---------|---------|
| `--gfp-dirname` | `gfp` | GFP channel folder name |
| `--cy-dirname` | `cy` | Cy channel folder name |
| `--roi-zip` | `roi.zip` | ROI zip filename, or full path for direct mode |
| `--channels` | both | Space-separated list of channel folders to measure |
| `-v / --verbose` | off | Verbose logging |

Works on macOS, Linux, and Windows. All dependencies are in the
repo-root `requirements.txt`.

---

### `golgi_signal_analysis.py` — Golgi localization quantification

Generates a Golgi mask from the Cy5 channel and measures how much GFP
signal falls inside vs. outside the Golgi, per Z slice, per cell.

For each cell, the script:
1. Computes an IsoData threshold on the Cy5 signal within the cell ROI
2. Applies that threshold to create a binary Golgi mask
3. Measures GFP intensity inside the Golgi mask and in the rest of the
   cell
4. Writes per-slice metrics to a CSV
5. Saves a binary mask TIFF and an ImageJ ROI zip for visual verification
   in FIJI

**Input:** any folder containing `gfp/` + `cy/` channel subfolders and
a `roi.zip` (one ROI per cell). Works with `Cropped/` or
`Background_Subtracted/` output.

**Output** (written to `<sample>/Results/`):

| File | Description |
|------|-------------|
| `golgi_signal_analysis.csv` | Per-slice metrics for all cells |
| `masks/<cell>_golgi_mask.tif` | Binary Golgi mask stack (Z, H, W) uint8 |
| `masks/<cell>_golgi_rois.zip` | ImageJ ROI zip — one polygon per Z slice |

**CSV columns:**

| Column | Description |
|--------|-------------|
| `cell` | Cell identifier (GFP filename stem, e.g. `gfp01`) |
| `slice` | Z-slice index (1-based) |
| `cy5_threshold` | IsoData threshold used to define the Golgi mask |
| `area_cell` | Total cell area in pixels |
| `area_inside` | Pixels inside the Golgi mask |
| `area_outside` | Pixels outside the Golgi mask (within the cell) |
| `gfp_mean_inside` | Mean GFP intensity inside the Golgi |
| `gfp_mean_outside` | Mean GFP intensity outside the Golgi |
| `gfp_integrated_inside` | Sum of GFP pixel values inside the Golgi |
| `gfp_integrated_outside` | Sum of GFP pixel values outside the Golgi |
| `gfp_integrated_total` | Sum of GFP pixel values across the whole cell |
| `fraction_gfp_inside` | `integrated_inside / integrated_total` |

Note: `area_inside + area_outside == area_cell` for every row. The mask
is always derived from the Cy5 channel regardless of the output filename.

**Usage:**

**Auto-discovery mode** — point at a root folder, the script finds all
samples:

```bash
python golgi_signal_analysis.py --input-dir /path/to/main_folder
```

**Direct mode** — when the `roi.zip` lives in a different folder from
the channel TIFFs (e.g. `roi.zip` in `Cropped/` but images in
`Background_Subtracted/`):

```bash
python golgi_signal_analysis.py \
    --input-dir /path/to/Background_Subtracted \
    --roi-zip /path/to/Cropped/roi.zip
```

**Verifying masks in FIJI:**

1. Open a Cy5 stack in FIJI
2. Go to `Analyze → Tools → ROI Manager → Open`, select
   `<cell>_golgi_rois.zip`
3. Scroll through Z — each slice should have a polygon outlining the
   Golgi signal
4. Alternatively, open `<cell>_golgi_mask.tif` as an overlay directly
   on top of the Cy5/GFP images

**Common options:**

| Flag | Default | Purpose |
|------|---------|---------|
| `--gfp-dirname` | `gfp` | GFP channel folder name |
| `--cy-dirname` | `cy` | Cy channel folder name |
| `--roi-zip` | `roi.zip` | ROI zip filename, or full path for direct mode |
| `-v / --verbose` | off | Verbose logging |

### `remove_body.py` — Erase a drawn ROI from paired GFP/Cy stacks

Zeroes out the pixels inside each image's ROI (e.g. a cell body outline),
on every Z slice of both channels — useful for isolating signal *outside*
a region, such as background or neurite signal once the cell body is
removed.

**Input:** a parent folder containing `gfp/`, `cy/`, and `roi/`
subfolders, where `roi/` has one ImageJ ROI `.zip` per image (e.g.
`01.zip`, `02.zip`, ...). Files across all three folders are paired
positionally after sorting by the numeric part of their filename — this
matches the numbering already used by this project's ROI-drawing output
(e.g. `roi_original/01.zip` aligned 1:1 with sequentially-renumbered
channel stacks).

**Output** (written to `<parent_dir>/Body_Removed/`):

```
Body_Removed/
    gfp/    gfp01.tif  gfp02.tif  ...
    cy/     cy01.tif   cy02.tif   ...
```

Erased pixels are set to `0`; everything outside the ROI is left
untouched. If a ROI zip contains more than one polygon, all of them are
unioned into a single erased region.

**Usage:**

```bash
python remove_body.py --input-dir /path/to/parent_dir
```

**Common options:**

| Flag | Default | Purpose |
|------|---------|---------|
| `--gfp-dirname` | `gfp` | GFP channel folder name |
| `--cy-dirname` | `cy` | Cy channel folder name |
| `--roi-dirname` | `roi` | ROI folder name (one `.zip` per image) |
| `--output-dirname` | `Body_Removed` | Output folder name, written under `--input-dir` |
| `-v / --verbose` | off | Verbose logging |

Works on macOS, Linux, and Windows. All dependencies are in the
repo-root `requirements.txt`.
