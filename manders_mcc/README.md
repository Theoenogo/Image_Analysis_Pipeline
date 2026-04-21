# Manders-MCC

Colocalization analysis tools for fluorescence microscopy images.  Computes
per-slice auto-thresholds, Pearson's correlation coefficient, and thresholded
Manders' colocalization coefficients (M1, M2) for paired GFP / CY image stacks.

## Scripts

| Script | Runs in | Speed | Dependencies |
|--------|---------|-------|--------------|
| `imagej_combined_analysis.py` | ImageJ / Fiji | Moderate | None (uses built-in ImageJ API) |
| `standalone_analysis.py` | Terminal (Python 3) | Fast | `numpy`, `tifffile` |

Both scripts produce identical output columns and use the same algorithms
(IsoData threshold, standard Pearson's, thresholded Manders').

## Expected folder structure

```
Background_Subtracted/
  gfp/
    gfp1.tif
    gfp2.tif
    ...
  cy/
    cy1.tif
    cy2.tif
    ...
```

Files are paired by number: `gfp1.tif` matches `cy1.tif`, etc.

## Usage

### ImageJ (combined script)

1. Open ImageJ / Fiji.
2. **Plugins > Macros > Run...** and select `imagej_combined_analysis.py`.
3. Choose your `Background_Subtracted` folder when prompted.
4. Results are written to `coloc_results.csv` and `thresholds_multislice.csv`
   inside that folder.

### Terminal (standalone script -- recommended for large datasets)

macOS / Linux:
```bash
pip install numpy tifffile
python standalone_analysis.py /path/to/Background_Subtracted
```

Windows:
```powershell
pip install numpy tifffile
python standalone_analysis.py C:\path\to\Background_Subtracted
```

Both scripts are pure Python and work identically on macOS, Linux,
and Windows.

## Output files

### `coloc_results.csv`

| gfp_image | cy_image | slice | gfp_threshold | cy_threshold | pearson | mandersA | mandersB |
|-----------|----------|-------|---------------|--------------|---------|----------|----------|

### `thresholds_multislice.csv`

Same format as the original thresholding script output, for backward
compatibility.

## Performance comparison

| Approach | Per-slice overhead | Window management | Plugin dependency |
|----------|-------------------|-------------------|-------------------|
| Original (threshold script + JACoP script) | ~2-3 s (sleeps + GUI) | Opens / leaks windows | JACoP |
| `imagej_combined_analysis.py` | ~0.5-2 s (pixel math in Jython) | No windows opened | None |
| `standalone_analysis.py` (NumPy) | ~0.01-0.05 s (vectorized) | N/A | numpy, tifffile |

## Methods

- **Thresholding** -- ImageJ's Default method (IsoData / iterative intermeans).
- **Pearson's r** -- Standard correlation coefficient over all pixel pairs in
  each slice.
- **Manders' M1 / M2** -- Thresholded colocalization coefficients.  M1 is the
  fraction of channel-A signal that colocalizes with channel B (both channels
  above their respective thresholds); M2 is the converse.
