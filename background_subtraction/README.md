# Background Subtraction

Per-cell, per-slice background subtraction. Generalization of the
ImageJ
[`Background_Subtraction_Tailored.ijm`](./imagej_reference/Background_Subtraction_Tailored.ijm)
macro: subtracts a tailored constant from each cropped single-cell
stack so only puncta above the cell's diffuse background remain —
input for `manders_mcc/` colocalization.

The original macro measured the mean inside the ROI **only on slice 0**
and reused that constant for the whole Z stack. This Python port
measures **per slice**, so Z-dependent haze is subtracted correctly.
For single-slice images the two behaviors are identical.

## What it does

For every `Cropped/` folder produced by `roi_cropping/`:

1. Walk the `gfp/` and `cy/` subfolders, sorting files numerically by
   the first digit run in the filename (`gfp1.tif < gfp2.tif < gfp10.tif`).
2. Open the combined `roi.zip`. Pair the *i*-th ROI with the *i*-th
   image pair (this matches `ROI_Select_Duplicate_TIFF_Loop.ijm`'s
   numbering, and is what `roi_cropping/` writes by construction).
3. For each image pair (per channel independently):
   - For every Z slice, compute the **mean intensity inside the ROI on
     that slice**.
   - Per slice, `subtract_value = mean * multiplier`, clamped to
     `[100, 5000]`.
   - Subtract each slice's value from that slice. Negative
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

Works on macOS, Linux, and Windows without platform-specific setup.

## Usage

```bash
python bg_subtract.py --input-dir /path/to/main_folder
```

The walker finds every `Cropped/` folder anywhere under `--input-dir`.

### Direct mode — one specific folder

If you want to background-subtract a single folder that isn't part of
the `roi_cropping/` walk (e.g. output from `extras/remove_body.py`),
use `--direct`. It treats `--input-dir` as the sample folder itself
instead of searching for `Cropped/` subfolders:

```bash
python bg_subtract.py --input-dir /path/to/some_folder --direct
```

`some_folder` must directly contain `gfp/`, `cy/`, and one of:

- a single `roi.zip` (one ROI per image, matched by index) — the
  original convention, or
- a `roi/` folder with one ROI file per image (`.roi` or `.zip`,
  matched positionally by the numeric part of the filename) — e.g.
  `extras/remove_body.py`'s `roi/` input. If a per-image file contains
  more than one ROI, they're unioned into a single mask.

`roi.zip` takes precedence if both are present. Output goes to a
sibling `Background_Subtracted/` folder next to `some_folder`, same as
normal mode.

The `roi/`-folder convention is also recognized automatically during
the regular recursive `Cropped/` walk (no `--direct` needed) — it's a
second accepted ROI source, not a replacement for `roi.zip`.

### Common options

| Flag | Default | Purpose |
|------|---------|---------|
| `--direct` | off | Treat `--input-dir` as the sample folder itself instead of walking for `Cropped/` subfolders. |
| `--roi-zip` | `roi.zip` | ROI zip filename to look for inside each sample. |
| `--roi-dirname` | `roi` | Fallback ROI folder name (one `.roi`/`.zip` per image) if `--roi-zip` isn't found. |
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

**Caveat:** the in-ROI mean is only a good stand-in for "diffuse
background" when the ROI doesn't contain a strong, spatially-concentrated
bright region. If it does (e.g. a bright Golgi/TGN signal in a
proinsulin or insulin channel), that region pulls the mean up, inflating
the subtraction constant — which then gets applied to the *whole* image,
disproportionately erasing real, dimmer signal far from that region
(puncta out in cell processes / near the plasma membrane). See
`bg_subtract_percentile.py` below for a fix.

## Alternative: percentile-based subtraction (`bg_subtract_percentile.py`)

Same algorithm shape — one flat constant per slice, clamped, applied to
the whole image — but computed from a **low percentile of the ROI**
(default: 10th) instead of the mean. A low percentile reflects the
typical dim background pixel regardless of how bright a small region
elsewhere in the ROI gets, since that region only ever occupies a
minority of the ROI's pixels — so it doesn't have the Golgi-inflation
problem above. Deliberately the smallest change that fixes it, not a new
algorithm: same inputs, same output layout, same `--gfp-multiplier` /
`--cy-multiplier` / `--floor` / `--ceiling` defaults as `bg_subtract.py`.

Writes to a *different* output folder
(`Background_Subtracted_Percentile/` by default) so it can sit alongside
`bg_subtract.py`'s own output for direct comparison — running it never
touches or overwrites the original method's results.

```bash
python bg_subtract_percentile.py --input-dir /path/to/main_folder --percentile 10
```

Accepts the same `--direct`, `--roi-zip`, `--roi-dirname` options as
`bg_subtract.py` (see [Common options](#common-options)), plus:

| Flag | Default | Purpose |
|------|---------|---------|
| `--percentile` | `10` | ROI percentile used in place of the mean. Lower = more conservative (subtracts less); higher moves back toward mean-like behavior and re-exposes the original problem as it approaches ~50. |
| `--output-dirname` | `Background_Subtracted_Percentile` | Output folder name, written as a sibling of `Cropped/`. |

Validated on real data (INS-1 endogenous proinsulin/insulin + CGA
datasets): with the defaults above, Golgi detection in a downstream
spot-counting step was essentially unchanged (~1% area difference vs.
the mean-based method), while the fraction of peripheral pixels zeroed
by the subtraction dropped sharply at distance from the Golgi (100–150px
band: 94% → 41%; 150–250px band: 92% → 54%) — i.e. real distal signal
that the mean-based method was erasing is preserved. `--percentile 10`
is a reasonable starting point, not a rigorously optimized value — tune
it for your own images the same way you'd tune any other threshold in
this pipeline: compare a few cells against the mean-based output at
matched display scaling before trusting a full batch run.

(An earlier spatially-local rolling-ball approach was tried and dropped:
it correctly decouples a bright Golgi from the rest of the image, but a
rolling ball can't remove a genuinely flat, uniform background — there's
no local trend for it to roll under — so it left far-field noise almost
entirely unsuppressed. The percentile fix above solves the actual
problem without that failure mode.)

## Parameter provenance

| Parameter | ImageJ source | This port |
|-----------|---------------|-----------|
| GFP multiplier | macro `getNumber(... 1.25)` L21 | `1.25` (same) |
| Cy multiplier | macro `getNumber(... 1.25)` L22 | `1.25` (same) |
| Lower clamp | macro `if (subtractValue < 100)` L132 | `100` (same) |
| Upper clamp | macro `if (subtractValue > 5000)` L133 | `5000` (same) |
| Measurement scope | macro `Measure` on active slice (L125) — slice 0 only | **per slice** (intentional generalization) |
| Subtraction scope | macro `for (s=1..nSlices)` (L138-142) — same constant for every slice | **per-slice value applied per slice** |

The only intentional deviation from the macro is the per-slice
measurement, which reduces to the macro's behavior on single-slice
inputs and is strictly more accurate on Z stacks.
