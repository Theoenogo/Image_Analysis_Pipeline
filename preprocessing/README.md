# Preprocessing

Python port of the MATLAB preprocessing pipeline that used to run before
`roi_drawing/`. Takes a folder of raw microscope TIFF slices, stacks them,
corrects GFP channel registration with a fixed XY offset, and runs
Richardson-Lucy deconvolution against measured PSFs.

Supports both **GFP/Cy** and **GFP/RFP** experiments. A sample normally
pairs GFP with one or the other, not both — pass `--cy-psf` and/or
`--rfp-psf` depending on which second channel this dataset uses. The GFP
XY offset is calibrated per pairing (see [GFP XY offset](#gfp-xy-offset)
below), so set `--gfp-offset` to match whichever channel you're running.

The original MATLAB scripts are preserved for reference in
[`matlab_reference/`](./matlab_reference).

## What it does

In a single command, this replaces MATLAB scripts A → F:

1. **Cleanup** — renames scope-output folders like `gfp1.abc123` or
   `gfp2.2026-04-10-16-34-59` to just `gfp1` / `gfp2`, and recursively
   deletes `*.scanprotocol` sidecar files.
2. **Stack + GFP XY offset** — for every `gfp*` folder under each sample,
   reads all single-page TIFFs, applies a fixed pixel shift to correct
   chromatic registration (default `(+5, −2)`, calibrated for GFP/Cy5 —
   see [GFP XY offset](#gfp-xy-offset)), and writes a multi-page uint16
   TIFF into `<sample>/Decon/gfp/<name>.tif`.
3. **Stack Cy** — same as step 2 for `cy*` folders, but with no offset.
4. **Stack RFP** — same as step 3 for `rfp*` folders, if present.
5. **Consolidate** — moves each `<sample>/Decon/` up to a single
   `<input-dir>/Decon/<sample>/` so every sample's stacks live under one
   tree.
6. **Richardson-Lucy deconvolution** — 30 iterations (matching the MATLAB
   pipeline's `DL2.RL(..., 30)`) against user-supplied PSF TIFFs, saved
   into `<input-dir>/Decon/Deconvoluted/<sample>/{gfp,cy,rfp}/<name>_decon.tif`.
   Only channels with a `--*-psf` supplied are deconvolved — GFP is
   always required, and at least one of `--cy-psf` / `--rfp-psf`.

Stacking always attempts all three channel prefixes; a sample with no
`rfp*`/`cy*` folders just produces zero stacks for that channel, so
GFP/Cy-only and GFP/RFP-only datasets both work without extra flags.

The `Deconvoluted/` tree is what `roi_drawing/` reads from next.

## GFP XY offset

The chromatic registration shift between GFP and the second channel
depends on which filter cube pairing was used:

| Pairing | `--gfp-offset` |
|---|---|
| GFP / Cy5 (default) | `5 -2` |
| GFP / RFP | `-3 1` |

Pass the value matching your dataset explicitly — it is not inferred from
which PSFs you supply:

```bash
python preprocess.py --input-dir ... --gfp-psf ... --rfp-psf ... \
    --gfp-offset -3 1
```

## EVOS microscope images

If your images were acquired on an EVOS microscope, both channels are
saved together in a single timestamped folder rather than in separate
`gfp{N}/` and `cy{N}/` subfolders. Run `extras/split_channels.py`
first to reorganize them into the layout below before running this
script.

EVOS images are saved as OME-TIFFs. The stacking step handles these
correctly — it reads only the first frame from any file that `tifffile`
reports as multi-frame (a common side-effect of OME metadata declaring
extra frames that were never written).

## Input layout

The simple 2-level layout (mirrors the MATLAB pipeline):

```
<input-dir>/
    sampleA/
        gfp1/                   ← individual .tif slices from the scope
        cy1/                    (GFP/Cy5 experiments)
        rfp1/                   (GFP/RFP experiments)
    sampleB/
        gfp2/
        cy2/
    ...
```

A sample normally has `cy*` *or* `rfp*`, not both. Both are stacked
opportunistically regardless (harmless if the folder isn't there — you
just get zero stacks for that channel), and only the channel(s) you pass
a `--*-psf` for are deconvolved.

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
            rfp/rfp1.tif        ← stacked, no shift (if present)
        sampleB/
            gfp/gfp2.tif
            cy/cy2.tif
        Deconvoluted/
            sampleA/
                gfp/gfp1_decon.tif  ← Richardson-Lucy output, uint16
                cy/cy1_decon.tif
                rfp/rfp1_decon.tif
            sampleB/
                gfp/gfp2_decon.tif
                cy/cy2_decon.tif
```

## Install

All dependencies come from the repo-root `requirements.txt`.

The default `dl2` engine matches the original MATLAB pipeline
bit-for-bit and requires:

- **A JDK (Java 8+)**:
  - macOS: `brew install openjdk`
  - Ubuntu/Debian: `sudo apt install default-jdk`
  - Windows: download from <https://adoptium.net/> (Temurin JDK),
    run the `.msi` installer, and make sure "Set JAVA_HOME" is
    checked during install.
- **Fiji** with the **DeconvolutionLab2** plugin installed:
  1. Download Fiji for your platform: <https://imagej.net/software/fiji/downloads>
  2. Open Fiji → `Help → Update... → Manage Update Sites`
  3. Enable **DeconvolutionLab2**, then click `Apply and Close`
  4. Click `Apply changes` in the updater.
  5. On Windows, extract the Fiji `.zip` to a short path like
     `C:\Fiji.app` (long paths with spaces can cause issues with
     DL2's internal path handling).
- **A display**: DL2 internally constructs Swing/AWT components even
  for programmatic calls, so it cannot run truly headless. Works on
  any local desktop (macOS, Linux, Windows). Will NOT work over
  SSH-without-X-forwarding or in a headless server/CI environment.

If you're running headless (server, SSH-only) or just don't want a
Fiji/JDK install, use the pure-Python fallback — `--engine scipy` needs
nothing beyond the standard requirements.

For GPU-accelerated deconvolution via PyTorch (MPS on Apple Silicon,
CUDA on NVIDIA), use `--engine torch` after:

```bash
pip install torch
```

Point the CLI at your Fiji installation:

macOS:
```bash
--fiji-dir /Applications/Fiji.app
# or: export FIJI_DIR=/Applications/Fiji.app
```

Windows:
```powershell
--fiji-dir C:\Fiji.app
# or: $env:FIJI_DIR = "C:\Fiji.app"
```

Linux:
```bash
--fiji-dir ~/Fiji.app
# or: export FIJI_DIR=~/Fiji.app
```

## Usage

### GFP / Cy5

macOS / Linux:
```bash
python preprocess.py \
    --input-dir /path/to/main_folder \
    --gfp-psf   /path/to/gfp_psf.tif \
    --cy-psf    /path/to/cy_psf.tif \
    --fiji-dir  /Applications/Fiji.app
```

Windows:
```powershell
python preprocess.py `
    --input-dir C:\path\to\main_folder `
    --gfp-psf   C:\path\to\gfp_psf.tif `
    --cy-psf    C:\path\to\cy_psf.tif `
    --fiji-dir  C:\Fiji.app
```

### GFP / RFP

Same as above with `--rfp-psf` instead of `--cy-psf`, and the RFP-specific
GFP offset (see [GFP XY offset](#gfp-xy-offset)):

```bash
python preprocess.py \
    --input-dir /path/to/main_folder \
    --gfp-psf   /path/to/gfp_psf.tif \
    --rfp-psf   /path/to/rfp_psf.tif \
    --gfp-offset -3 1 \
    --fiji-dir  /Applications/Fiji.app
```

Channels deconvolve **concurrently** by default — the script launches one
worker thread per channel sharing a single Fiji JVM, cutting wall-clock
time roughly in half (or a third, for all three channels) compared to
running them sequentially. Pass `--sequential-channels` to disable.

Independently, pass `--jobs N` to also deconvolve multiple images
*within* each channel concurrently — the same concurrency model as
[`deconvolve_folder.py --jobs`](#standalone-single-folder-deconvolution).
The two combine multiplicatively: `--jobs 3` on a GFP/RFP run (2 channels)
can launch up to 2 × 3 = 6 Fiji processes at once. Lower `--jobs` or add
`--sequential-channels` on memory-constrained machines.

```bash
python preprocess.py --input-dir ... --gfp-psf ... --cy-psf ... --jobs 3
```

For `--engine dl2`, each concurrent job — whether from `--jobs` or from
channel-level parallelism — launches its own independent Fiji subprocess.
DL2's multi-page TIFF writer isn't atomic, and an early version of this
pipeline had a race condition where quitting Fiji as soon as the output
file *appeared* could truncate the later slices of the stack, especially
under concurrent load. `deconvolve_dl2.py`'s `deconvolve_file()` now
guards against this for every call, concurrent or not: unique per-call
image titles, `-port0` to stop concurrent launches collapsing onto one
shared JVM, an in-macro slice-count re-verification loop before Fiji is
allowed to quit, and a Python-side page-count check that raises a clear
`RuntimeError` if a result still comes out short. Any deconvolution
failure — including a truncated file caught by that check — aborts the
whole run instead of silently producing a short batch.

### Common options

- `--iterations 30` — Richardson-Lucy iteration count. Default `30`
  matches the MATLAB pipeline.
- `--gfp-offset X Y` — XY pixel shift applied to every GFP slice. Default
  `5 -2` (GFP/Cy5); use `-3 1` for GFP/RFP. Pass `0 0` to disable. See
  [GFP XY offset](#gfp-xy-offset).
- `--engine dl2|scipy|torch` — deconvolution engine.
  - `dl2` (default) calls the DeconvolutionLab2 Java plugin via PyImageJ.
    Matches the original MATLAB pipeline bit-for-bit. Needs Fiji + JDK.
  - `scipy` is a pure-Python FFT-based Richardson-Lucy on CPU — no Fiji,
    no JVM. Not numerically identical to DL2 but close.
  - `torch` is the same pure-Python RL, but in PyTorch for
    MPS / CUDA acceleration.
- `--fiji-dir /path/to/Fiji.app` — only used with `--engine dl2`.
  Also reads `$FIJI_DIR` from the environment.
- `--torch-device cpu|mps|cuda` — only used with `--engine torch`;
  default `mps`.
- `--no-crop-psf` — disable automatic PSF cropping. By default the PSF
  is cropped to its signal bounding box (+ 30 px margin) before
  deconvolution. This is safe and dramatically speeds up large PSFs
  captured at full sensor resolution (e.g. 1328×1048 → ~60×60).
- `--psf-crop-margin N` — pixel margin around the PSF signal region
  when auto-cropping (default `30`). Increase if you want to preserve
  more of the outer diffraction rings.
- `--sequential-channels` — run channels back-to-back instead of
  concurrently. Useful for debugging or on very memory-constrained
  machines.
- `--jobs N` / `-j N` — deconvolve N images within each channel
  concurrently (default `1`). See the note above on how this combines
  with channel-level concurrency.
- `--skip-cleanup`, `--skip-stacking`, `--skip-consolidate`,
  `--skip-deconvolution` — run subsets of the pipeline. Useful for
  re-running just the deconvolution after updating a PSF.

### Apple Silicon (M1/M2/M3/M4)

- `--engine dl2` works fine under Apple Silicon; make sure you're using
  a `arm64` JDK (e.g. `brew install openjdk`).
- The `scipy` fallback goes through Apple's Accelerate framework via
  pocketfft and is fast in practice.
- For GPU-accelerated Python RL via Metal:

```bash
python preprocess.py --input-dir ... --gfp-psf ... --cy-psf ... \
    --engine torch --torch-device mps
```

### Windows

- Install a 64-bit JDK from <https://adoptium.net/> (Temurin). Check
  "Set JAVA_HOME" during install. Verify with `java -version` in
  PowerShell.
- Download the Windows `.zip` of Fiji from <https://imagej.net/software/fiji/downloads>,
  extract to a short path like `C:\Fiji.app`. The launcher is
  `ImageJ-win64.exe` inside that directory.
- Install DL2 in Fiji the same way as on Mac (Help → Update →
  Manage Update Sites → DeconvolutionLab2).
- Use `--fiji-dir C:\Fiji.app` or set `$env:FIJI_DIR = "C:\Fiji.app"`
  in your PowerShell profile.
- The `scipy` engine needs no extra setup and is the easiest way to
  get started on Windows.
- For GPU-accelerated Python RL via CUDA:

```powershell
pip install torch
python preprocess.py --input-dir ... --gfp-psf ... --cy-psf ... `
    --engine torch --torch-device cuda
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

## Standalone single-folder deconvolution

If you just want to deconvolve one folder of TIFF stacks against one PSF —
no stacking, no consolidation, no gfp/cy pairing — use
`deconvolve_folder.py` instead of the full pipeline:

```bash
python deconvolve_folder.py \
    --input-dir  /path/to/some_folder \
    --psf        /path/to/psf.tif \
    --output-dir /path/to/output_folder \
    --fiji-dir   /Applications/Fiji.app
```

Every `*.tif` directly inside `--input-dir` (non-recursive) is deconvolved
and written to `<output-dir>/<name>_decon.tif`. If `--output-dir` is
omitted, outputs go to `<input-dir>_decon` next to the input folder.

It accepts the same `--iterations`, `--engine`, `--fiji-dir`,
`--torch-device`, `--no-crop-psf`, and `--psf-crop-margin` options as
`preprocess.py` (see [Common options](#common-options) above) — it's the
same deconvolution code, just pointed at a single flat folder instead of
the pipeline's `<sample>/<channel_group>/` layout.

Pass `--jobs N` (default `1`) to deconvolve multiple images concurrently
instead of one at a time:

```bash
python deconvolve_folder.py \
    --input-dir  /path/to/some_folder \
    --psf        /path/to/psf.tif \
    --output-dir /path/to/output_folder \
    --fiji-dir   /Applications/Fiji.app \
    --jobs 3
```

For `--engine dl2`, each job launches its own Fiji process (JVM + GUI
window), so `--jobs` genuinely runs several deconvolutions side by side —
2-4 is a reasonable starting point depending on CPU cores and RAM. For
`--engine scipy`/`torch` it uses threads instead, which mainly helps if
the computation isn't already saturating all CPU cores (or the GPU, for
`torch`).

## Algorithm details

Richardson-Lucy implemented as standard multiplicative updates:

```
estimate_{n+1} = estimate_n * conv(image / conv(estimate_n, psf), mirror(psf))
```

- **`dl2`** engine: calls DeconvolutionLab2 directly via PyImageJ, with
  the same arguments the MATLAB pipeline used:
  `DL2.RL(image, psf, 30.0, '-out stack short')`. Image I/O is handled
  by ImageJ itself (`IJ.openImage` / `IJ.saveAsTiff`), so the on-disk
  format matches ImageJ's convention exactly.
- **`scipy` / `torch`** engines: FFT-based Python RL using
  `scipy.signal.fftconvolve` or `torch.fft.fftn`. Both produce identical
  results to each other up to floating-point noise. Output is rescaled
  to the uint16 range only when RL produces values outside
  `[0, 65535]`, matching the MATLAB behavior. PSFs are normalized to
  sum to 1 before use.

## Parameter provenance

The hardcoded defaults mirror the MATLAB pipeline:

| Parameter              | MATLAB source                                   | Default |
|------------------------|-------------------------------------------------|---------|
| GFP XY offset (vs Cy5) | `C2_green_stack_offset_Cy5_recursive.m` L43–44  | `5, -2` |
| GFP XY offset (vs RFP) | user-supplied, not from MATLAB reference        | `-3, 1` |
| Cy XY offset           | `D2_cy5_stack_no_offset_recursive.m`            | none    |
| RFP XY offset          | none (same treatment as Cy)                     | none    |
| R-L iterations         | `F_Cy5_runDeconvolutionLab2_parallel.m` L91     | `30`    |
| Output dtype           | `F_*` scripts, `uint16()` cast                  | uint16  |
| Uncompressed TIFF      | `F_*` scripts, `'Compression', 'none'`          | yes     |
| Deconvolution engine   | `DL2.RL(image, psf, 30, '-out stack short')`    | same (default) |
