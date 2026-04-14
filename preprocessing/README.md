# Preprocessing

Python port of the MATLAB preprocessing pipeline that used to run before
`roi_drawing/`. Takes a folder of raw microscope TIFF slices, stacks them,
corrects GFP/Cy channel registration with a fixed XY offset, and runs
Richardson-Lucy deconvolution against measured PSFs.

The original MATLAB scripts are preserved for reference in
[`matlab_reference/`](./matlab_reference).

## What it does

In a single command, this replaces MATLAB scripts A → F:

1. **Cleanup** — renames scope-output folders like `gfp1.abc123` or
   `gfp2.2026-04-10-16-34-59` to just `gfp1` / `gfp2`, and recursively
   deletes `*.scanprotocol` sidecar files.
2. **Stack + GFP XY offset** — for every `gfp*` folder under each sample,
   reads all single-page TIFFs, applies a fixed `(+5, −2)` pixel shift to
   correct chromatic registration, and writes a multi-page uint16 TIFF
   into `<sample>/Decon/gfp/<name>.tif`.
3. **Stack Cy** — same as step 2 for `cy*` folders, but with no offset.
4. **Consolidate** — moves each `<sample>/Decon/` up to a single
   `<input-dir>/Decon/<sample>/` so every sample's stacks live under one
   tree.
5. **Richardson-Lucy deconvolution** — 30 iterations (matching the MATLAB
   pipeline's `DL2.RL(..., 30)`) against user-supplied PSF TIFFs, saved
   into `<input-dir>/Decon/Deconvoluted/<sample>/{gfp,cy}/<name>_decon.tif`.

The `Deconvoluted/` tree is what `roi_drawing/` reads from next.

## Input layout

The simple 2-level layout (mirrors the MATLAB pipeline):

```
<input-dir>/
    sampleA/
        gfp1/                   ← individual .tif slices from the scope
        cy1/
        rfp1/                   (ignored by default; add a third stacking
                                 pass if you need RFP)
    sampleB/
        gfp2/
        cy2/
    ...
```

Folders can be named with or without the post-dot suffix (`gfp1`,
`gfp1.abc123`, or `gfp2.2026-04-10-16-34-59`) — the cleanup step
normalizes them all to just `gfp1`, `gfp2`, etc.

### Nested / deeper layouts

Any depth of nesting is supported. All stages walk the tree
recursively and find sample folders wherever they live. For example:

```
<input-dir>/
    experiment1/
        replicate1/
            sampleA/
                gfp1/
                cy1/
            sampleB/
                gfp2/
                cy2/
        replicate2/
            sampleC/
                gfp3/
                cy3/
    experiment2/
        sampleD/
            gfp4/
            cy4/
```

When there are intermediate group/replicate folders, the consolidate
step flattens the sample path into the consolidated folder name using
``__`` as the separator. E.g. the sample at
``<input-dir>/experiment1/replicate1/sampleA`` becomes
``<input-dir>/Decon/experiment1__replicate1__sampleA/``. For the flat
2-level layout the consolidated name is just the sample folder name
(``sampleA``), unchanged from the MATLAB behavior.

Because the deconvolution step writes per-sample output
(``Deconvoluted/<sample>/{gfp,cy}/...``), channel folder names do
*not* need to be unique across the whole experiment — two samples
can each have a ``gfp1`` folder without collision.

## Output layout

After a full run:

```
<input-dir>/
    Decon/
        sampleA/
            gfp/gfp1.tif        ← stacked, XY-shifted, pre-deconvolution
            cy/cy1.tif          ← stacked, no shift, pre-deconvolution
        sampleB/
            gfp/gfp2.tif
            cy/cy2.tif
        Deconvoluted/
            sampleA/
                gfp/gfp1_decon.tif  ← Richardson-Lucy output, uint16
                cy/cy1_decon.tif
            sampleB/
                gfp/gfp2_decon.tif
                cy/cy2_decon.tif
```

## Install

All dependencies come from the repo-root `requirements.txt`. If you want
the optional torch/MPS deconvolution backend on Apple Silicon, install
PyTorch as well:

```bash
pip install torch
```

## Usage

```bash
python preprocess.py \
    --input-dir /path/to/main_folder \
    --gfp-psf   /path/to/gfp_psf.tif \
    --cy-psf    /path/to/cy_psf.tif
```

### Common options

- `--iterations 30` — Richardson-Lucy iteration count. Default `30`
  matches the MATLAB pipeline.
- `--gfp-offset X Y` — XY pixel shift applied to every GFP slice. Default
  `5 -2`. Pass `0 0` to disable.
- `--backend scipy|torch` — FFT backend for deconvolution. `scipy`
  (default) is CPU-based and fast on Apple Silicon via Accelerate.
  `torch` enables MPS / CUDA.
- `--torch-device cpu|mps|cuda` — only used with `--backend torch`;
  default `mps`.
- `--skip-cleanup`, `--skip-stacking`, `--skip-consolidate`,
  `--skip-deconvolution` — run subsets of the pipeline. Useful for
  re-running just the deconvolution after updating a PSF.

### Apple Silicon (M1/M2/M3/M4)

The default `scipy` backend is typically fast enough — SciPy's FFT goes
through Apple's Accelerate framework, which is well-tuned for these
machines. If you want to try GPU-accelerated deconvolution via Metal:

```bash
python preprocess.py --input-dir ... --gfp-psf ... --cy-psf ... \
    --backend torch --torch-device mps
```

### Example: re-deconvolve only

If you've already stacked and consolidated and just want to redo
deconvolution with a different PSF or iteration count:

```bash
python preprocess.py \
    --input-dir /path/to/main_folder \
    --gfp-psf   /path/to/new_gfp_psf.tif \
    --cy-psf    /path/to/new_cy_psf.tif \
    --iterations 50 \
    --skip-cleanup --skip-stacking --skip-consolidate
```

## Algorithm details

Richardson-Lucy implemented as standard multiplicative updates:

```
estimate_{n+1} = estimate_n * conv(image / conv(estimate_n, psf), mirror(psf))
```

Both backends use FFT-based convolution (`scipy.signal.fftconvolve` or
`torch.fft.fftn`) and produce identical results up to floating-point
noise. Output is rescaled to the uint16 range only when Richardson-Lucy
produces values outside `[0, 65535]`, matching the MATLAB behavior.

PSFs are normalized to sum to 1 before use.

## Parameter provenance

The hardcoded defaults mirror the MATLAB pipeline:

| Parameter              | MATLAB source                                   | Default |
|------------------------|-------------------------------------------------|---------|
| GFP XY offset          | `C2_green_stack_offset_Cy5_recursive.m` L43–44  | `5, -2` |
| Cy XY offset           | `D2_cy5_stack_no_offset_recursive.m`            | none    |
| R-L iterations         | `F_Cy5_runDeconvolutionLab2_parallel.m` L91     | `30`    |
| Output dtype           | `F_*` scripts, `uint16()` cast                  | uint16  |
| Uncompressed TIFF      | `F_*` scripts, `'Compression', 'none'`          | yes     |
