"""Richardson-Lucy deconvolution — port of MATLAB F_GFP / F_Cy5 scripts.

The MATLAB pipeline used DeconvolutionLab2 with ``DL2.RL(image, psf, 30)`` —
classic Richardson-Lucy with 30 iterations against a user-supplied measured
PSF. This module supports three interchangeable engines:

- **dl2** (default): calls the actual DeconvolutionLab2 Java plugin via
  PyImageJ, so results match the MATLAB pipeline bit-for-bit. Requires
  a local Fiji install with DL2 plugin + a JDK. See ``deconvolve_dl2.py``.
- **scipy**: ``scipy.signal.fftconvolve`` on CPU. Pure Python — no JVM,
  no Fiji, no plugin. Numerically close to DL2 but not identical
  (boundary handling / convergence differ).
- **torch**: same algorithm as scipy but in PyTorch, which enables MPS
  acceleration on Apple Silicon or CUDA on NVIDIA.

The scipy and torch engines produce identical output up to floating-point
noise. DL2 is the reference and should be preferred whenever it is
available.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional

import numpy as np
import tifffile
from scipy.signal import fftconvolve

log = logging.getLogger(__name__)


Engine = Literal["dl2", "scipy", "torch"]
TorchDevice = Literal["cpu", "mps", "cuda"]


@dataclass
class DeconvolutionConfig:
    num_iter: int = 30
    engine: Engine = "dl2"
    torch_device: TorchDevice = "mps"
    fiji_dir: Optional[Path] = None  # only used when engine == "dl2"
    crop_psf: bool = True
    psf_crop_margin: int = 30


# ---------------------------------------------------------------------------
# PSF auto-crop
# ---------------------------------------------------------------------------


def _autocrop_psf(
    psf: np.ndarray,
    margin: int = 30,
    threshold_frac: float = 0.001,
) -> np.ndarray:
    """Crop a PSF to its signal bounding box plus *margin* pixels per side.

    Measured PSFs are often captured at full sensor resolution but contain
    signal in a small central region. Cropping removes the vast black
    border so FFT-based deconvolution operates on a much smaller volume.

    *threshold_frac* is the fraction of the PSF max used to distinguish
    signal from background. 0.1% is conservative enough to preserve faint
    outer diffraction rings while ignoring camera dark-current noise.
    """
    peak = float(psf.max())
    if peak <= 0:
        log.warning("PSF max is <= 0; skipping auto-crop")
        return psf

    threshold = peak * threshold_frac
    nonzero = np.argwhere(psf > threshold)
    if nonzero.size == 0:
        log.warning("No PSF pixels above %.3f%% of max; skipping auto-crop",
                     threshold_frac * 100)
        return psf

    mins = nonzero.min(axis=0)
    maxs = nonzero.max(axis=0)

    slices = []
    for dim in range(psf.ndim):
        lo = max(0, int(mins[dim]) - margin)
        hi = min(psf.shape[dim], int(maxs[dim]) + margin + 1)
        slices.append(slice(lo, hi))

    cropped = psf[tuple(slices)].copy()
    old_shape = "x".join(str(s) for s in psf.shape)
    new_shape = "x".join(str(s) for s in cropped.shape)
    log.info("Auto-cropped PSF from %s to %s (margin=%d, threshold=%.1f%% of max)",
             old_shape, new_shape, margin, threshold_frac * 100)
    return cropped


def _save_cropped_psf(psf: np.ndarray, suffix: str = "") -> Path:
    """Write a cropped PSF to a temp file and return its path."""
    import tempfile
    fd, path = tempfile.mkstemp(suffix=f"_psf{suffix}.tif")
    os.close(fd)
    tifffile.imwrite(path, psf.astype(np.uint16), imagej=True,
                     photometric="minisblack")
    return Path(path)


# ---------------------------------------------------------------------------
# Core Richardson-Lucy
# ---------------------------------------------------------------------------


def _richardson_lucy_scipy(
    image: np.ndarray,
    psf: np.ndarray,
    num_iter: int,
) -> np.ndarray:
    """CPU Richardson-Lucy using scipy's FFT-based convolution.

    Works on 2D or 3D arrays (image and psf must share dimensionality).
    Returns a float32 array the same shape as ``image``.
    """
    image = image.astype(np.float32, copy=False)
    psf = psf.astype(np.float32, copy=True)
    psf_sum = psf.sum()
    if psf_sum <= 0:
        raise ValueError("PSF sum is non-positive; cannot normalize.")
    psf /= psf_sum
    psf_mirror = psf[tuple(slice(None, None, -1) for _ in psf.shape)]

    estimate = np.clip(image, a_min=1e-8, a_max=None).copy()
    eps = 1e-8
    for _ in range(num_iter):
        blurred = fftconvolve(estimate, psf, mode="same")
        np.maximum(blurred, eps, out=blurred)
        ratio = image / blurred
        correction = fftconvolve(ratio, psf_mirror, mode="same")
        estimate *= correction
        np.maximum(estimate, 0.0, out=estimate)
    return estimate


def _richardson_lucy_torch(
    image: np.ndarray,
    psf: np.ndarray,
    num_iter: int,
    device: TorchDevice,
) -> np.ndarray:
    """Torch Richardson-Lucy. FFT-based, uses the requested device."""
    try:
        import torch
    except ImportError as err:
        raise RuntimeError(
            "torch backend requested but PyTorch is not installed. "
            "Install with `pip install torch`, or use --backend scipy."
        ) from err

    if device == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError(
            "MPS requested but not available on this machine. "
            "Use --torch-device cpu or install a PyTorch build with MPS."
        )
    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but not available on this machine.")

    torch_device = torch.device(device)

    # Normalize PSF
    psf_np = psf.astype(np.float32)
    psf_np = psf_np / psf_np.sum()

    # Pad PSF to image shape, centered at origin with ifftshift semantics so
    # the FFT convolution matches a 'same' centered convolution.
    img_shape = image.shape
    pad_psf = np.zeros(img_shape, dtype=np.float32)
    slices_dst = tuple(slice(0, s) for s in psf_np.shape)
    pad_psf[slices_dst] = psf_np
    # Shift so PSF center is at index 0.
    shifts = tuple(-(s // 2) for s in psf_np.shape)
    pad_psf = np.roll(pad_psf, shifts, axis=tuple(range(pad_psf.ndim)))

    img_t = torch.as_tensor(image.astype(np.float32), device=torch_device)
    psf_t = torch.as_tensor(pad_psf, device=torch_device)

    otf = torch.fft.fftn(psf_t)
    otf_conj = otf.conj()

    estimate = torch.clamp(img_t, min=1e-8).clone()
    eps = 1e-8
    for _ in range(num_iter):
        blurred = torch.fft.ifftn(torch.fft.fftn(estimate) * otf).real
        blurred = torch.clamp(blurred, min=eps)
        ratio = img_t / blurred
        correction = torch.fft.ifftn(torch.fft.fftn(ratio) * otf_conj).real
        estimate = estimate * correction
        estimate = torch.clamp(estimate, min=0.0)
    return estimate.detach().to("cpu").numpy()


def richardson_lucy(
    image: np.ndarray,
    psf: np.ndarray,
    num_iter: int = 30,
    backend: Literal["scipy", "torch"] = "scipy",
    torch_device: TorchDevice = "cpu",
) -> np.ndarray:
    """Run pure-Python Richardson-Lucy deconvolution.

    This is the numpy/torch path used by ``engine="scipy"`` and
    ``engine="torch"``. The DL2 engine bypasses this function entirely
    and calls the Java plugin directly (see ``deconvolve_dl2.py``).

    Parameters
    ----------
    image, psf
        Arrays of matching dimensionality (2D or 3D).
    num_iter
        Iteration count. The original MATLAB pipeline uses 30.
    backend
        ``"scipy"`` (CPU) or ``"torch"`` (any torch device).
    torch_device
        Only used when ``backend="torch"``.
    """
    if image.ndim != psf.ndim:
        raise ValueError(
            f"image.ndim={image.ndim} does not match psf.ndim={psf.ndim}"
        )
    if backend == "scipy":
        return _richardson_lucy_scipy(image, psf, num_iter)
    if backend == "torch":
        return _richardson_lucy_torch(image, psf, num_iter, torch_device)
    raise ValueError(f"Unknown backend: {backend}")


# ---------------------------------------------------------------------------
# Stack I/O helpers
# ---------------------------------------------------------------------------


def _load_tiff_stack(path: Path) -> np.ndarray:
    """Load a multi-page TIFF as a 3D (Z, H, W) float32 array."""
    arr = tifffile.imread(str(path))
    if arr.ndim == 2:
        arr = arr[np.newaxis, ...]
    return arr.astype(np.float32)


def _save_uint16_stack(path: Path, stack: np.ndarray) -> None:
    """Rescale to uint16 (matching MATLAB) and save as an uncompressed multi-page TIFF.

    MATLAB behavior: if min < 0 or max > 65535, rescale to [0, 65535],
    otherwise clip to uint16 as-is.
    """
    min_val = float(stack.min())
    max_val = float(stack.max())
    if min_val < 0 or max_val > 65535:
        if max_val > min_val:
            rescaled = (stack - min_val) / (max_val - min_val) * 65535.0
        else:
            rescaled = np.zeros_like(stack)
    else:
        rescaled = stack
    out = np.clip(np.rint(rescaled), 0, 65535).astype(np.uint16)
    tifffile.imwrite(str(path), out, imagej=True, photometric="minisblack", compression=None)


# ---------------------------------------------------------------------------
# Per-channel orchestration
# ---------------------------------------------------------------------------


def _deconvolve_file_python(
    input_path: Path,
    output_path: Path,
    psf: np.ndarray,
    config: DeconvolutionConfig,
) -> None:
    """Per-file RL for the scipy / torch engines (pre-loaded PSF)."""
    image = _load_tiff_stack(input_path)
    if image.ndim != psf.ndim:
        raise ValueError(
            f"Dimension mismatch: image {image.shape} vs PSF {psf.shape}. "
            "Both must have the same number of axes."
        )
    backend = "scipy" if config.engine == "scipy" else "torch"
    result = richardson_lucy(
        image, psf,
        num_iter=config.num_iter,
        backend=backend,
        torch_device=config.torch_device,
    )
    _save_uint16_stack(output_path, result)


def _iter_channel_jobs(
    decon_folder: Path,
    channel_group: str,
) -> list[tuple[Path, Path]]:
    """Collect ``(input_tif, output_tif)`` pairs for every stack to deconvolve.

    Walks every ``<decon_folder>/<sample>/<channel_group>/*.tif`` pair and
    maps it to its matching
    ``<decon_folder>/Deconvoluted/<sample>/<channel_group>/<name>_decon.tif``.
    """
    deconvoluted_root = decon_folder / "Deconvoluted"
    jobs: list[tuple[Path, Path]] = []
    for channel_dir in sorted(decon_folder.rglob(channel_group)):
        if not channel_dir.is_dir():
            continue
        if deconvoluted_root in channel_dir.parents or channel_dir == deconvoluted_root:
            continue
        sample_rel = channel_dir.parent.relative_to(decon_folder)
        target_dir = deconvoluted_root / sample_rel / channel_group
        for stack_path in sorted(channel_dir.glob("*.tif")):
            if stack_path.name.endswith("_decon.tif"):
                continue
            out_path = target_dir / f"{stack_path.stem}_decon.tif"
            jobs.append((stack_path, out_path))
    return jobs


def deconvolve_channel(
    decon_folder: Path,
    psf_path: Path,
    channel_group: str,
    config: DeconvolutionConfig,
) -> list[Path]:
    """Deconvolve every stack in ``<decon_folder>/<sample>/<channel_group>/*.tif``.

    Mirrors the MATLAB F scripts: for each ``<channel_group>`` folder found
    recursively under ``decon_folder`` (typically
    ``<main>/Decon/<sample>/<group>/``), deconvolves every ``.tif`` stack
    inside it and writes the result to
    ``<decon_folder>/Deconvoluted/<sample>/<channel_group>/<name>_decon.tif``.

    The per-sample output layout means samples keep their channel stacks
    grouped together, and two samples that happen to share a channel folder
    name (e.g. both have ``gfp1``) don't overwrite each other.

    Parameters
    ----------
    decon_folder
        The consolidated ``Decon`` folder, typically ``<main_folder>/Decon``.
    psf_path
        Measured PSF TIFF (same dimensionality as the stacks, typically 3D).
    channel_group
        Folder name to match, e.g. ``"gfp"`` or ``"cy"``.
    config
        Algorithm parameters (iterations, engine, device).
    """
    decon_folder = Path(decon_folder)
    psf_path = Path(psf_path)
    if not decon_folder.is_dir():
        raise NotADirectoryError(decon_folder)
    if not psf_path.is_file():
        raise FileNotFoundError(psf_path)

    deconvoluted_root = decon_folder / "Deconvoluted"
    deconvoluted_root.mkdir(parents=True, exist_ok=True)

    jobs = _iter_channel_jobs(decon_folder, channel_group)
    if not jobs:
        log.info("No %s stacks to deconvolve under %s", channel_group, decon_folder)
        return []

    log.info("Found %d %s stack(s) to deconvolve (engine=%s)",
             len(jobs), channel_group, config.engine)

    outputs: list[Path] = []

    # --- PSF auto-crop (benefits all engines) ---
    effective_psf_path = psf_path
    cropped_psf_tmp: Path | None = None
    cropped_psf_arr: np.ndarray | None = None

    if config.crop_psf:
        raw_psf = _load_tiff_stack(psf_path)
        cropped = _autocrop_psf(raw_psf, margin=config.psf_crop_margin)
        if cropped.shape != raw_psf.shape:
            if config.engine == "dl2":
                cropped_psf_tmp = _save_cropped_psf(
                    cropped, suffix=f"_{channel_group}")
                effective_psf_path = cropped_psf_tmp
            else:
                cropped_psf_arr = cropped
        else:
            if config.engine != "dl2":
                cropped_psf_arr = raw_psf

    try:
        if config.engine == "dl2":
            from .deconvolve_dl2 import deconvolve_file as _dl2_deconvolve_file
            for input_path, output_path in jobs:
                output_path.parent.mkdir(parents=True, exist_ok=True)
                log.info("Deconvolving (DL2) %s", input_path)
                _dl2_deconvolve_file(
                    input_path, output_path, effective_psf_path,
                    num_iter=config.num_iter,
                    fiji_dir=config.fiji_dir,
                )
                outputs.append(output_path)
                log.info("Wrote %s", output_path)
        else:
            psf = cropped_psf_arr if cropped_psf_arr is not None else _load_tiff_stack(psf_path)
            log.info("Loaded PSF %s with shape %s", psf_path, psf.shape)
            for input_path, output_path in jobs:
                log.info("Deconvolving (%s) %s", config.engine, input_path)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                _deconvolve_file_python(input_path, output_path, psf, config)
                outputs.append(output_path)
                log.info("Wrote %s", output_path)
    finally:
        if cropped_psf_tmp is not None:
            try:
                cropped_psf_tmp.unlink()
            except OSError:
                pass

    return outputs
