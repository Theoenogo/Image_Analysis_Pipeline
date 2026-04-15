# Image Analysis Pipeline

A modular Python pipeline for analyzing widefield fluorescence microscopy images (.tif). Designed for quantifying colocalization between fluorescence channels in biological imaging data.

---

## Pipeline Overview

Tools should be run in the following order:

| Step | Folder | Description |
|------|--------|-------------|
| 1 | `preprocessing/` | Stack microscope TIFF slices, correct GFP XY offset, Richardson-Lucy deconvolution |
| 2 | `roi_drawing/` | Automated ROI detection around cells using Cellpose (writes to `roi_original/`) |
| 3 | *(manual)* ImageJ `Color_Merge_Automated_PreloadROIs_Adjust.ijm` | User refines ROIs and saves the edited set to `roi/` |
| 4 | `roi_cropping/` | Crop each ROI's single-cell stack from the decon images and re-origin the ROI coords |
| 5 | `background_subtraction/` | Per-cell background subtraction (mean-inside-ROI × multiplier, clamped) |
| 6 | `manders_mcc/` | Colocalization analysis (Manders' coefficients + Pearson's) |

---

## Requirements

- Python 3.12
- All dependencies listed in `requirements.txt`

---

## Setup

**1. Clone the repository**
```bash
git clone git@github.com:Theoenogo/Image_Analysis_Pipeline.git
cd Image_Analysis_Pipeline
```

**2. Create a virtual environment**
```bash
python3.12 -m venv venv
source venv/bin/activate
```

**3. Install dependencies**
```bash
python3.12 -m pip install -r requirements.txt
```

> Activate the virtual environment (`source venv/bin/activate`) at the start of every session before running any scripts.

---

## Tools

### Preprocessing (`preprocessing/`)

Stacks per-slice microscope TIFFs into multi-page stacks, applies a
fixed XY offset to the GFP channel to correct for channel registration,
and runs Richardson-Lucy deconvolution against user-supplied PSFs.
This is the Python port of the original MATLAB preprocessing pipeline.

**Run:**
```bash
cd preprocessing
python preprocess.py \
    --input-dir /path/to/main_folder \
    --gfp-psf /path/to/gfp_psf.tif \
    --cy-psf /path/to/cy_psf.tif
```

See `preprocessing/README.md` for full usage and options.

---

### ROI Drawing (`roi_drawing/`)

Automatically detects cells that are positive in both the GFP and Cy5 channels of a fluorescence image pair, draws ROIs around each cell's full membrane, and exports ImageJ-compatible ROI `.zip` files.

Built around [Cellpose](https://github.com/MouseLand/cellpose), which handles punctate cytoplasmic signal where classical threshold + watershed pipelines break down.

**Run:**
```bash
cd roi_drawing
python roi_detect.py
```

See `roi_drawing/README.md` for full usage and options.

---

### Manual ROI editing (ImageJ — kept as-is)

Between `roi_drawing/` and `roi_cropping/`, the user runs the ImageJ macro
`Color_Merge_Automated_PreloadROIs_Adjust.ijm` (preserved in
`roi_cropping/matlab_reference/`). It loads each image pair plus the
auto-generated ROIs from `roi_original/` so the user can add, remove,
or refine ROIs by hand, then saves the edited set to a sibling `roi/`
folder. This is the only manual step in the pipeline.

---

### ROI Cropping (`roi_cropping/`)

Reads the user-edited ROI zips and the decon stacks and produces one
cropped single-cell stack per ROI in flat `Cropped/{gfp,cy}/` folders,
with all pixels outside the ROI zeroed and the ROI re-origined to the
crop-local frame. Replaces the old MATLAB
`H_roi_duplicate_image_all_channels_Recursive.m` + ImageJ
`ROI_Select_Duplicate_TIFF_Loop.ijm` two-step flow with a single pass.

**Run:**
```bash
cd roi_cropping
python roi_crop.py --input-dir /path/to/main_folder
```

See `roi_cropping/README.md` for full usage and options.

---

### Background Subtraction (`background_subtraction/`)

For each cropped cell, measures the mean intensity inside the ROI on
slice 0, multiplies by a user-set multiplier (default 1.25, clamped to
`[100, 5000]`), and subtracts that constant from every Z slice. Writes
to `Background_Subtracted/{gfp,cy}/`, which is what `manders_mcc/`
reads. Python port of `Background_Subtraction_Tailored.ijm`.

**Run:**
```bash
cd background_subtraction
python bg_subtract.py --input-dir /path/to/main_folder
```

See `background_subtraction/README.md` for full usage and options.

---

### Manders MCC (`manders_mcc/`)

Computes per-slice auto-thresholds, Pearson's correlation coefficient, and thresholded Manders' colocalization coefficients (M1, M2) for paired GFP / CY image stacks.

Two scripts are available:
- `standalone_analysis.py` — run directly from the terminal
- `imagej_combined_analysis.py` — run inside Fiji/ImageJ

**Run (standalone):**
```bash
cd manders_mcc
python standalone_analysis.py
```

See `manders_mcc/README.md` for full usage and options.

---

## Batch runner (`run_experiment.py`)

A convenience wrapper that chains the three post-ImageJ-edit stages
(`roi_cropping → background_subtraction → manders_mcc`) across every
sample in an experiment. Use it when you've finished the manual ROI
editing in Fiji and want to go from `<sample>/{gfp,cy,roi}/` all the
way to final colocalization CSVs in one command.

**Run:**
```bash
python run_experiment.py --input-dir /path/to/main_folder
```

The individual stage scripts (`roi_crop.py`, `bg_subtract.py`,
`standalone_analysis.py`) are untouched and still usable on their own
when you want to iterate on a single stage. Stages can also be skipped
individually in the batch runner via `--skip-roi-cropping`,
`--skip-background-subtraction`, `--skip-manders-mcc`.

Preprocessing + ROI drawing are intentionally **not** bundled in here —
both need per-experiment inputs (PSFs, Cellpose parameters), and the
manual Fiji ROI edit sits between them, so chaining them with the
post-edit stages wouldn't be useful.

---

## Adding New Tools

When adding a new pipeline step:
1. Create a new subfolder (e.g. `preprocessing/`)
2. Add any new dependencies with `python3.12 -m pip install <package>`
3. Update the master requirements file: `python3.12 -m pip freeze > requirements.txt`
4. Update the Pipeline Overview table in this README

---

## Contributing

This is a private lab repository. Contact the repository owner before making changes.
