# Background Subtraction

Python port of the ImageJ
[`Background_Subtraction_Tailored.ijm`](./imagej_reference/Background_Subtraction_Tailored.ijm)
macro. Subtracts a per-cell constant from every cropped single-cell stack so
only puncta above the cell's diffuse background remain — input for
`manders_mcc/` colocalization.

## What it does

For every `Cropped/` folder produced by `roi_cropping/`:

1. Walk the `gfp/` and `cy/` subfolders, sorting files numerically by
   the first digit run in the filename (`gfp1.tif < gfp2.tif < gfp10.tif`).
2. Open the combined `roi.zip`. Pair the *i*-th ROI with the *i*-th
   image pair (this matches `ROI_Select_Duplicate_TIFF_Loop.ijm`'s
   numbering, and is what `roi_cropping/` writes by construction).
3. For each image pair (per channel independently):
   - Compute the **mean intensity inside the ROI on slice 0**
     (matches the ImageJ macro: it calls `Measure` once on the active
     slice, not across Z).
   - `subtract_value = mean * multiplier`, clamped to `[100, 5000]`.
   - Subtract that constant from every pixel of every Z slice. Negative
     intermediates are clipped at zero before the cast back to uint16.
4. Save into a sibling `Background_Subtracted/{gfp,cy}/` folder with
   the same filename as the source.

The output layout is exactly what `manders_mcc/` reads from next.

## Input layout

```
<sample>/Cropped/
    gfp/gfp1.tif, gfp2.tif, ...
    cy/cy1.tif, cy2.tif, ...
    roi.zip
```

## Output layout

```
<sample>/Background_Subtracted/
    gfp/gfp1.tif, gfp2.tif, ...
    cy/cy1.tif, cy2.tif, ...
```

## Install

All dependencies come from the repo-root `requirements.txt`.

## Usage

```bash
python bg_subtract.py --input-dir /path/to/main_folder
```

The walker finds every `Cropped/` folder anywhere under `--input-dir`.

### Common options

| Flag | Default | Purpose |
|------|---------|---------|
| `--gfp-multiplier` | `1.25` | Multiplier on the GFP per-cell mean. Matches the ImageJ macro default. |
| `--cy-multiplier` | `1.25` | Multiplier on the Cy per-cell mean. |
| `--floor` | `100` | Lower clamp on the subtraction value. |
| `--ceiling` | `5000` | Upper clamp on the subtraction value. |
| `-v / --verbose` | off | Verbose logging. |

## Why "mean inside ROI × 1.25"?

After `roi_cropping/` zeros out pixels outside the ROI, the in-ROI mean
is the cell's **diffuse cytoplasmic intensity**. Subtracting
`mean × 1.25` removes that diffuse background and leaves only puncta
that are at least 25% brighter than the cell average — which is what
`manders_mcc/` then thresholds and colocalizes.

The clamp (`[100, 5000]`) prevents pathological cases (very dim or very
bright cells) from producing extreme subtractions. Both the clamp and
the multiplier match the ImageJ macro defaults verbatim.

## Parameter provenance

| Parameter | ImageJ source | Default |
|-----------|---------------|---------|
| GFP multiplier | macro `getNumber(... 1.25)` L21 | `1.25` |
| Cy multiplier | macro `getNumber(... 1.25)` L22 | `1.25` |
| Lower clamp | macro `if (subtractValue < 100)` L132 | `100` |
| Upper clamp | macro `if (subtractValue > 5000)` L133 | `5000` |
| Measurement scope | macro `Measure` on active slice (L125) | slice 0 |
| Subtraction scope | macro `for (s=1..nSlices)` (L138-142) | every slice |
