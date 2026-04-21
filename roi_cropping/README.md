# ROI Cropping

Python port of the MATLAB ``H_roi_duplicate_image_all_channels_Recursive.m``
plus the ImageJ ``ROI_Select_Duplicate_TIFF_Loop.ijm`` two-step flow,
collapsed into a single pass.

The original MATLAB/ImageJ scripts are preserved for reference in
[`matlab_reference/`](./matlab_reference).

## What it does

Given a sample folder containing the per-FOV decon stacks plus the
user-edited ROI zips (one zip per FOV), produces one cropped single-cell
stack per ROI in flat per-channel folders ready for the next pipeline
stage.

In a single command, this replaces the duplicate-then-crop dance:

1. **Walk** the input directory for sample folders containing all of
   `gfp/`, `cy/`, and `roi/` subfolders.
2. **Match** each ROI zip to its image pair by sorted numeric index in
   the filename (`01.zip` ↔ `gfp1`/`cy1`, etc.).
3. **For every ROI in every zip**, compute the tight bounding box, crop
   both gfp and cy stacks to that bbox, and zero out pixels outside the
   ROI polygon (matches ImageJ's "Make Inverse → Clear stack").
4. **Re-origin** the ROI polygon coordinates to the bbox-local frame so
   the saved per-cell ROI overlays correctly on the cropped image.
5. **Save** sequentially-numbered crops + a combined `roi.zip` for
   downstream use.

The intermediate "duplicate each image N times" step from the MATLAB
flow was a MATLAB/ImageJ plumbing workaround and is **not needed in
Python** — we iterate over ROIs directly.

## Input layout

The input is the output of `roi_drawing/` after the user has manually
refined ROIs in the ImageJ
`Color_Merge_Automated_PreloadROIs_Adjust.ijm` macro. The edited ROIs
live in `roi/`, and the unedited automated set is preserved in
`roi_original/`:

```
<sample>/
    gfp/gfp1_decon.tif           ← from preprocessing/
    gfp/gfp2_decon.tif
    cy/cy1_decon.tif
    cy/cy2_decon.tif
    roi_original/01.zip          ← from roi_drawing/ (untouched)
    roi_original/02.zip
    roi/01.zip                   ← user-edited via Color_Merge ImageJ macro
    roi/02.zip
```

Sample discovery is recursive — any folder anywhere under `--input-dir`
that contains all three of `gfp/`, `cy/`, and `roi/` is treated as a
sample. This handles the nested
`Decon/Deconvoluted/<sample>/` layout from `preprocessing/` automatically.

## Output layout

```
<sample>/Cropped/
    gfp/
        gfp1.tif        ← single-cell stack, ROI-shaped, all-zero outside
        gfp2.tif
        ...
    cy/
        cy1.tif
        cy2.tif
        ...
    roi/
        1.roi           ← per-cell ROI in bbox-local coords
        2.roi
        ...
    roi.zip             ← combined zip in image-paired order (single ROI per cell)
```

Numbering is sequential across all FOVs in the sample: 4 ROIs in pair
01 + 3 ROIs in pair 02 yields `gfp1.tif … gfp7.tif`. The combined
`roi.zip` enumerates ROIs in the same order, so downstream stages can
pair `gfpN.tif`/`cyN.tif` with the Nth ROI directly.

## Install

All dependencies come from the repo-root `requirements.txt`. The
relevant packages are `tifffile`, `roifile`, `scikit-image`, `numpy`.

Works on macOS, Linux, and Windows without platform-specific setup.

## Usage

```bash
python roi_crop.py --input-dir /path/to/main_folder
```

If your channel folders are named differently (e.g. `gfp1`, `cy5`):

```bash
python roi_crop.py \
    --input-dir /path/to/main_folder \
    --gfp-dirname gfp1 \
    --cy-dirname cy5
```

### Common options

| Flag | Default | Purpose |
|------|---------|---------|
| `--gfp-dirname` | `gfp` | Name of the GFP channel folder inside each sample. |
| `--cy-dirname` | `cy` | Name of the Cy channel folder. |
| `--roi-dirname` | `roi` | Name of the user-edited ROI folder. |
| `--gfp-prefix` | `gfp` | Filename prefix for cropped GFP outputs. |
| `--cy-prefix` | `cy` | Filename prefix for cropped Cy outputs. |
| `-v / --verbose` | off | Verbose logging. |

## Algorithm details

For each ROI:

1. **Bounding box** — `(xmin, ymin)` to `(xmax, ymax)` over the ROI
   polygon, clipped to the image extents.
2. **Crop** — slice both stacks to the bbox over all Z planes.
3. **Mask** — rasterize the polygon into a bbox-local boolean mask via
   `skimage.draw.polygon`. Pixels outside the polygon get zeroed in
   every Z plane.
4. **Re-origin ROI** — subtract `(xmin, ymin)` from every polygon
   vertex, rebuild as an `ImagejRoi` so it overlays correctly on the
   cropped image.
5. **Write** — preserve source dtype (uint16 from decon), save
   uncompressed multi-page TIFF.

Output dtype matches the source decon stacks (typically uint16),
identical to the ImageJ macro's behavior.
