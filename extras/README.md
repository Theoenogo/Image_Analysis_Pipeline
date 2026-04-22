# Extras

Standalone analysis scripts that are useful but not part of the main
pipeline (preprocessing → ROI drawing → cropping → background
subtraction → colocalization). Each script is self-contained and can
be run independently.

## Scripts

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

```bash
python measure_roi_signal.py --input-dir /path/to/main_folder
```

To measure only one channel:

```bash
python measure_roi_signal.py --input-dir /path/to/main_folder --channels cy
```

To point at a specific subfolder (e.g. `Cropped/` directly):

```bash
python measure_roi_signal.py --input-dir /path/to/sample/Cropped
```

**Common options:**

| Flag | Default | Purpose |
|------|---------|---------|
| `--gfp-dirname` | `gfp` | GFP channel folder name |
| `--cy-dirname` | `cy` | Cy channel folder name |
| `--roi-zip` | `roi.zip` | ROI zip filename |
| `--channels` | both | Space-separated list of channel folders to measure |
| `-v / --verbose` | off | Verbose logging |

Works on macOS, Linux, and Windows. All dependencies are in the
repo-root `requirements.txt`.
