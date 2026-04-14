# Image Analysis Pipeline

A modular Python pipeline for analyzing widefield fluorescence microscopy images (.tif). Designed for quantifying colocalization between fluorescence channels in biological imaging data.

---

## Pipeline Overview

Tools should be run in the following order:

| Step | Folder | Description |
|------|--------|-------------|
| 1 | *(coming soon)* | Pre-processing / image preparation |
| 2 | `roi_drawing/` | Automated ROI detection around cells using Cellpose |
| 3 | *(coming soon)* | Intermediate processing |
| 4 | `manders_mcc/` | Colocalization analysis (Manders' coefficients + Pearson's) |

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

## Adding New Tools

When adding a new pipeline step:
1. Create a new subfolder (e.g. `preprocessing/`)
2. Add any new dependencies with `python3.12 -m pip install <package>`
3. Update the master requirements file: `python3.12 -m pip freeze > requirements.txt`
4. Update the Pipeline Overview table in this README

---

## Contributing

This is a private lab repository. Contact the repository owner before making changes.
