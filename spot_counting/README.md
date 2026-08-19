# Spot Counting

Counts granule puncta in one channel and scores whether each one overlaps a
punctum in the paired channel. Automates the manual protocol: pick up to 20
puncta per cell in the granule-marker channel (proinsulin or mature insulin),
look at the paired channel (e.g. CGA), call overlap yes/no, report a percent
overlap per cell.

Different from `manders_mcc/`: that measures intensity *correlation* over whole
images (Pearson's r, Manders' M1/M2). This one works on **discrete spots** and
answers "does this particular granule have CGA on it?"

## The Golgi problem

The strong perinuclear signal is proinsulin at the Golgi. It is not punctate and
must never be counted as a granule. It is excluded automatically:

1. IsoData-threshold the reference channel inside the cell.
2. Label connected components. Anything at or above `--golgi-min-area` (default
   300 px) is Golgi — granule puncta are a handful of pixels, the Golgi ribbon is
   hundreds to thousands, so area alone separates them.
3. Dilate by `--golgi-dilate` px so the bright rim can't spawn false puncta.

The same exclusion is applied to the CGA channel, which also has TGN signal.

**Every cell gets a QC overlay** (`Results/qc/<cell>_spots.png`) showing the
excluded region and the scored spots. Look at these before trusting a batch.

## Marker rules — you must set `--marker`

`--marker` is **required and never inferred from the image**. Both markers show
perinuclear signal plus surrounding puncta; the difference is what you *intend*
to count, not something measurable in the pixels. Run the script once per marker.

| `--marker` | Rule | Why |
|---|---|---|
| `proinsulin` | Golgi blob excluded; juxta-Golgi puncta **kept** | Immature granules bud at the Golgi and move outward, so the puncta right next to the ribbon are real objects |
| `insulin` | Golgi excluded **and** puncta within `--min-golgi-distance` (default 15 px of the Golgi edge) dropped | Mature granules are the ones out in the processes/arms and along the plasma membrane |

The marker and every parameter it implies are written into each output row, so
no number is ever separated from the rule that produced it.

If proinsulin puncta are being picked too far out into the processes for your
liking, add `--max-golgi-distance N` — the opposite of `--min-golgi-distance`:
it drops any punctum farther than N px from the Golgi instead of closer than
N px, restricting the marker to clearly near-Golgi/immature puncta. Not tied
to a `--marker` preset (no validated default exists yet) — tune it by eye:
the QC overlay draws the cutoff as a dashed orange ring, so you can see
directly which puncta it would include or exclude before committing to a
value.

## How a cell is analyzed

1. Max-intensity-project both channels over Z (analysis is 2D).
2. Percentile-stretch each projection to `[0, 1]` using percentiles measured
   **inside the cell mask**. This is what lets one `--spot-threshold` work
   across images instead of needing retuning per image.
3. Build the Golgi exclusion mask (above).
4. Detect puncta with Laplacian-of-Gaussian (`skimage.feature.blob_log`) in both
   channels, keeping centres inside the cell and outside the Golgi. Detection
   runs on the unmodified image and centres are filtered afterwards — blanking
   the Golgi first would create a hard edge the LoG filter answers with spurious
   responses.
5. If more than `--spots-per-cell` valid spots remain, draw a **random** sample
   using `--seed`. Random rather than brightest-first: ranking by intensity
   would bias toward large or bright granules and inflate the overlap estimate.
6. For each counted spot, find the nearest spot in the other channel. Overlap =
   yes if that distance is within `--match-radius`, which is what absorbs XY
   drift between channels.

## Input layout

Any stage that has `gfp/`, `cy/` and `roi.zip` — normally `Background_Subtracted/`,
since removing diffuse haze makes puncta stand out:

```
<sample>/
    gfp/gfp01.tif  gfp02.tif  …    (one multi-page stack per cell)
    cy/cy01.tif    cy02.tif   …    (matched 1-to-1 with gfp/)
    roi.zip                         (ImageJ ROI zip, one ROI per cell)
```

Files are paired positionally after sorting by the numeric part of the filename,
and the *i*-th ROI goes with the *i*-th image pair — the same convention the rest
of the pipeline uses.

## Usage

```bash
python spot_count.py --input-dir /path/to/proinsulin_data --marker proinsulin
python spot_count.py --input-dir /path/to/insulin_data    --marker insulin
```

The walker finds every qualifying sample folder under `--input-dir`. Use
`--direct` to analyze one specific folder without walking.

**Verify which physical channel (`gfp`/`cy`) actually holds the granule
marker before running — don't assume.** `gfp`/`cy` are generic microscope
filter-channel names, not protein labels; which antibody went on which
fluorophore is a fact about your staining protocol, and it can differ
between imaging sessions even for the "same" experiment type. Getting
this backwards doesn't error — it silently runs the whole analysis with
the two channels' roles swapped (Golgi exclusion, chance-level baseline,
everything). Default is `--reference-channel gfp`; add
`--reference-channel cy` if the granule marker is in Cy instead.

If the ROI zip lives outside `--input-dir` — e.g. pointing at
`Background_Subtracted/` while `roi.zip` is still in the sibling
`Cropped/` folder that produced it (`background_subtraction` never copies
it over) — pass the full path via `--roi-zip` together with `--direct`:

```bash
python spot_count.py --input-dir /path/to/Background_Subtracted \
    --roi-zip /path/to/Cropped/roi.zip --marker proinsulin --direct
```

## Options

| Flag | Default | Notes |
|---|---|---|
| `--marker` | **required** | `proinsulin` or `insulin` |
| `--reference-channel` | `gfp` | Channel holding the granule marker. **Set deliberately.** |
| `--spots-per-cell` | `20` | Matches the prior manual protocol |
| `--seed` | `0` | Makes the random draw reproducible |
| `--spot-min-sigma` / `--spot-max-sigma` | `1.0` / `3.0` | Puncta size range in px |
| `--spot-threshold` | `0.05` | LoG cutoff. Lower catches dimmer granules, higher rejects noise |
| `--golgi-min-area` | `300` | px; area at or above which a bright region is Golgi |
| `--golgi-dilate` | `3` | px halo around the exclusion zone |
| `--match-radius` | `3` | px; XY-drift tolerance for calling overlap |
| `--min-golgi-distance` | per `--marker` | Override the marker preset; `0` disables |
| `--max-golgi-distance` | `0` (disabled) | Drop puncta farther than this from the Golgi. Not set by `--marker` presets — tune by eye |
| `--chance-permutations` | `200` | Draws per cell for the chance-level baseline (see below). Cheap; raise for a tighter estimate |
| `--pixel-size-nm` | none | If given, the three distance flags are read in nm |
| `--gfp-dirname` / `--cy-dirname` | `gfp` / `cy` | Channel folder names |
| `--roi-zip` | `roi.zip` | ROI zip filename inside each sample, or a full path elsewhere (requires `--direct`) |
| `--direct` | off | Treat `--input-dir` as one sample |
| `-v` | off | Verbose logging |

## Output

Written to `<sample>/Results/`:

| File | Contents |
|---|---|
| `spot_counts.csv` | One row per cell: spots detected/counted/overlapping, percent overlap, chance-level baseline, Golgi area, thresholds, and every parameter used |
| `spot_details.csv` | One row per counted spot: position, sigma, both channels' intensities, nearest-neighbour distance, overlap yes/no |
| `qc/<cell>_spots.png` | Overlay — yellow = Golgi excluded, green = overlap, red = no overlap, blue `+` = other-channel spots. Title includes the chance-level baseline |
| `qc/<cell>_spots.zip` | ImageJ ROIs of the counted spots; names end in `_ov`/`_no` |

Plus `spot_summary.csv` at the root of `--input-dir`: one row per sample with
mean percent overlap, SEM, mean overlap-above-chance, SEM, and n cells — the
table to run the proinsulin-vs-insulin comparison on. **Use the
above-chance columns for that comparison**, not raw percent overlap — see
below.

`spot_details.csv` records the other channel's intensity at every counted spot,
so overlap can be re-thresholded post-hoc without re-running the detection.

## Chance-level correction

A dense "other" channel can produce a high `percent_overlap` from pure
spatial coincidence, with zero true colocalization — e.g. if the other
channel's spots average 9 px apart, roughly 29% of random locations will
land within a 3 px `--match-radius` of one just by crowding. This is not
hypothetical: on real data, the same match radius gave a ~3% chance
baseline in one dataset and a ~29% baseline in another, purely because
one channel had far denser detections than the other — datasets acquired
on different days can differ this much in staining intensity/density.
Comparing raw `percent_overlap` between cells or datasets with different
"other"-channel density is therefore misleading on its own.

For every cell, the script also computes what overlap rate you'd expect
from chance alone: it takes the actual counted reference spots, then
repeatedly scatters the same *number* of "other"-channel spots uniformly
at random within the same valid region (same cell shape, same Golgi
exclusion) and recomputes overlap each time (`--chance-permutations`
draws, default 200). The columns this produces:

| Column | Meaning |
|---|---|
| `other_spots_detected` | How many spots were detected in the other channel — the density driving the baseline |
| `expected_overlap_chance_pct` | Mean overlap rate from `--chance-permutations` random draws — the "coincidence floor" |
| `chance_std_pct` | Standard deviation across those draws |
| `overlap_above_chance_pct` | `percent_overlap - expected_overlap_chance_pct` — the part of the observed overlap not explained by density alone |

**Use `overlap_above_chance_pct` (and its summary-level mean/SEM) for
comparisons across cells or datasets** — it's the number that isolates
real colocalization signal from how crowded each channel happens to be.

## Tuning and validation

Parameter tuning matters more here than code correctness. In order:

0. **Confirm which channel is which.** Before anything else, verify with
   whoever ran the staining/imaging which of `gfp`/`cy` actually holds the
   granule marker for *this* dataset, and pass `--reference-channel`
   accordingly. Don't assume it matches a different dataset from the same
   project — imaging sessions on different dates can differ.
1. **Eyeball one sample.** Run with `--direct`, open a `qc/*_spots.png`, and
   check: the Golgi is inside the yellow contour, circles sit on real puncta,
   and nothing is circled on noise or on the Golgi rim.
2. **Tune.** Golgi not fully excluded → raise `--golgi-min-area` or
   `--golgi-dilate`. Dim granules missed → lower `--spot-threshold`. Noise
   counted → raise it. Circles visibly the wrong size → adjust the sigma range.
   Counted spots landing too far out in the processes for a near-Golgi marker
   → add `--max-golgi-distance N`, adjusting N against the dashed orange ring
   until it matches where you'd stop calling puncta "immature."
3. **Check against manual counts.** Drag `qc/<cell>_spots.zip` onto the image in
   Fiji for a few cells scored by hand and compare spot-for-spot. If percent
   overlap is systematically off, the detection parameters are wrong.
4. **Confirm the CGA control.** CGA sorts from the TGN into the budding
   granule and stays with it, so CGA/proinsulin and CGA/insulin
   `overlap_above_chance_pct` (not raw `percent_overlap` — see
   [Chance-level correction](#chance-level-correction)) should be
   statistically indistinguishable **when each marker is scored under its
   own default rule** (proinsulin unrestricted, insulin peripheral-only).

   Don't "fix" the comparison by forcing both arms to `--min-golgi-distance
   0` and treating a flipped result as proof of an artifact — for an
   insulin antibody that cross-reacts with unprocessed proinsulin (common;
   confirm against your reagent's spec), that same-rule test isn't neutral.
   Near the Golgi, "insulin channel" signal in that setup is largely
   cross-reactive proinsulin, not mature insulin, so it will tend to
   resemble the proinsulin arm's near-Golgi puncta (including in CGA
   overlap) — a flip there reflects antibody specificity, not a bug in the
   peripheral restriction. Verified on a real 20230214/20230327 dataset
   pair: default rules gave p=0.89 (no difference, control passes); forcing
   both to `--min-golgi-distance 0` gave p=0.033 (insulin significantly
   higher) — expected once you account for cross-reactivity, not evidence
   the default rule is wrong. If your insulin antibody is confirmed
   proinsulin-specific-free (doesn't cross-react), this caveat doesn't
   apply and the same-rule test is a legitimate check.
5. **Reproducibility.** Re-run with the same `--seed` and confirm identical output.

Only then batch the full dataset.

## Install

All dependencies come from the repo-root `requirements.txt`. Works on macOS,
Linux and Windows without platform-specific setup.
