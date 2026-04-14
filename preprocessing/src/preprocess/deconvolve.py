"""Richardson-Lucy deconvolution — port of MATLAB F_GFP / F_Cy5 scripts.

The MATLAB pipeline used DeconvolutionLab2 with ``DL2.RL(image, psf, 30)`` —
classic Richardson-Lucy with 30 iterations against a user-supplied measured
PSF. This module implements the same algorithm in Python with two
interchangeable FFT backends:

- **scipy** (default): ``scipy.signal.fftconvolve`` on CPU. Uses pocketfft,
  which on Apple Silicon is routed through the Accelerate framework and
  is fast enough in practice for typical widefield stacks.
- **torch** (optional): same algorithm in PyTorch, which enables MPS
  acceleration on Apple Silicon ``Mac`` or CUDA on NVIDIA. Loaded lazily
  so the scipy path has no torch dependency.

Both backends produce the same output up to floating-point noise.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import tifffile
from scipy.signal import fftconvolve

log = logging.getLogger(__name__)


Backend = Literal["scipy", "torch"]
TorchDevice = Literal["cpu", "mps", "cuda"]


@dataclass
class DeconvolutionConfig:
    num_iter: int = 30
    backend: Backend = "scipy"
    torch_device: TorchDevice = "cpu"


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
    backend: Backend = "scipy",
    torch_device: TorchDevice = "cpu",
) -> np.ndarray:
    """Run Richardson-Lucy deconvolution.

    Parameters
    ----------
    image, psf
        Arrays of matching dimensionality (2D or 3D).
    num_iter
        Iteration count. The original MATLAB pipeline uses 30.
    backend
        ``"scipy"`` (CPU, default) or ``"torch"`` (any torch device).
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
    tifffile.imwrite(str(path), out, photometric="minisblack", compression=None)


# ---------------------------------------------------------------------------
# Per-channel orchestration
# ---------------------------------------------------------------------------


def deconvolve_channel(
    decon_folder: Path,
    psf_path: Path,
    channel_group: str,
    config: DeconvolutionConfig,
) -> list[Path]:
    """Deconvolve every stack in ``<decon_folder>/**/<channel_group>/*.tif``.

    Mirrors the MATLAB F scripts: for each ``<channel_group>`` folder found
    recursively under ``decon_folder`` (typically ``<main>/Decon/<sample>/<group>/``),
    deconvolves every ``.tif`` stack inside it, saves the result as
    ``<name>_decon.tif`` in the same folder, then moves all ``*_decon.tif``
    to ``<decon_folder>/Deconvoluted/<channel_group>/``.

    Parameters
    ----------
    decon_folder
        The consolidated ``Decon`` folder, typically ``<main_folder>/Decon``.
    psf_path
        Measured PSF TIFF (same dimensionality as the stacks, typically 3D).
    channel_group
        Folder name to match, e.g. ``"gfp"`` or ``"cy"``.
    config
        Algorithm parameters (iterations, backend, device).
    """
    decon_folder = Path(decon_folder)
    psf_path = Path(psf_path)
    if not decon_folder.is_dir():
        raise NotADirectoryError(decon_folder)
    if not psf_path.is_file():
        raise FileNotFoundError(psf_path)

    psf = _load_tiff_stack(psf_path)
    log.info("Loaded PSF %s with shape %s", psf_path, psf.shape)

    target_dir = decon_folder / "Deconvoluted" / channel_group
    target_dir.mkdir(parents=True, exist_ok=True)

    # First pass: collect every stack that would be deconvolved and check for
    # flat-output filename collisions across samples (which can happen when
    # nested sample folders share channel-folder names like "gfp1").
    stacks_to_process: list[Path] = []
    seen_output_names: dict[str, Path] = {}
    collisions: list[tuple[Path, Path]] = []
    for channel_dir in sorted(decon_folder.rglob(channel_group)):
        if not channel_dir.is_dir():
            continue
        if (decon_folder / "Deconvoluted") in channel_dir.parents:
            continue
        for stack_path in sorted(channel_dir.glob("*.tif")):
            if stack_path.name.endswith("_decon.tif"):
                continue
            output_name = f"{stack_path.stem}_decon.tif"
            if output_name in seen_output_names:
                collisions.append((seen_output_names[output_name], stack_path))
            else:
                seen_output_names[output_name] = stack_path
            stacks_to_process.append(stack_path)

    if collisions:
        collision_detail = "\n".join(
            f"  - {a}\n    collides with {b}" for a, b in collisions
        )
        raise ValueError(
            f"Flat output layout would overwrite {len(collisions)} "
            f"file(s) in Deconvoluted/{channel_group}/ because multiple "
            "samples produced stacks with the same name. "
            "Rename channel folders so each has a unique name across the "
            f"whole experiment (e.g. gfp1, gfp2, ...), or keep your data "
            "in a flat 2-level layout.\n\nCollisions:\n"
            f"{collision_detail}"
        )

    outputs: list[Path] = []
    for stack_path in stacks_to_process:
        log.info("Deconvolving %s", stack_path)
        image = _load_tiff_stack(stack_path)
        if image.ndim != psf.ndim:
            # Bail loudly rather than silently misbehaving.
            raise ValueError(
                f"Dimension mismatch: image {image.shape} vs PSF {psf.shape}. "
                "Both must have the same number of axes."
            )
        result = richardson_lucy(
            image, psf,
            num_iter=config.num_iter,
            backend=config.backend,
            torch_device=config.torch_device,
        )
        output_name = f"{stack_path.stem}_decon.tif"
        final_path = target_dir / output_name
        _save_uint16_stack(final_path, result)
        outputs.append(final_path)
        log.info("Wrote %s", final_path)

    return outputs
