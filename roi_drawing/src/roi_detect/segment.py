"""Cellpose wrapper: preprocess the pan-cell channel and run segmentation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class SegmentationParams:
    model_name: str = "cyto3"
    diameter: float | None = None  # None means let Cellpose auto-estimate
    blur_sigma: float = 2.0
    flow_threshold: float = 0.4
    cellprob_threshold: float = 0.0
    # Number of dynamics iterations. 0 means "Cellpose default" (~200).
    # Increase to 2000+ for cells with long thin processes so the flow
    # field has time to propagate all the way down the processes.
    niter: int = 0
    # Minimum detected object size in pixels. Drops tiny fragments.
    min_size: int = 15
    use_gpu: bool = True


def irregular_shapes_preset(params: "SegmentationParams") -> "SegmentationParams":
    """Tune parameters for cells with non-round shapes or long processes.

    Good starting point for INS-1 832/13 cells, neurons, or any
    adherent culture where the cells have arms / protrusions.
    """
    params.blur_sigma = 1.0
    params.flow_threshold = 0.6
    params.cellprob_threshold = -3.0
    params.niter = 4000
    return params


def preprocess_composite(
    channel_a: np.ndarray, channel_b: np.ndarray, blur_sigma: float = 1.0
) -> np.ndarray:
    """Create a composite segmentation input from two channels.

    Each channel is independently percentile-normalized and blurred, then
    the element-wise maximum is taken. This ensures the segmentation mask
    covers the full extent of signal in *either* channel — so ROIs won't
    clip GFP processes that extend beyond the Cy5 cell boundary.
    """
    a_norm = preprocess_for_cellpose(channel_a, blur_sigma=blur_sigma)
    b_norm = preprocess_for_cellpose(channel_b, blur_sigma=blur_sigma)
    return np.maximum(a_norm, b_norm)


def preprocess_for_cellpose(image: np.ndarray, blur_sigma: float = 2.0) -> np.ndarray:
    """Prepare a pan-cell channel for Cellpose.

    1. Percentile-normalize to [0, 1] using the 1st and 99.5th percentiles
       (robust contrast stretch, immune to hot pixels).
    2. Apply Gaussian blur with the given sigma to bridge punctate signal
       into continuous cell bodies.

    Returns a float32 array with the same shape as the input.
    """
    from skimage.filters import gaussian

    image = image.astype(np.float32)
    lo, hi = np.percentile(image, (1.0, 99.5))
    if hi <= lo:
        normalized = np.zeros_like(image, dtype=np.float32)
    else:
        normalized = np.clip((image - lo) / (hi - lo), 0.0, 1.0)

    if blur_sigma > 0:
        normalized = gaussian(normalized, sigma=blur_sigma, preserve_range=True)

    return normalized.astype(np.float32)


def _cellpose_major_version() -> int:
    """Return the major version of the installed Cellpose package."""
    try:
        import cellpose  # type: ignore

        return int(cellpose.__version__.split(".")[0])
    except Exception:
        return 3  # assume v3 if we can't determine


class CellposeSegmenter:
    """Lazy Cellpose model wrapper.

    The Cellpose model is expensive to construct (loads weights from disk),
    so build one instance and reuse it across every image in a batch.
    """

    def __init__(self, params: SegmentationParams):
        self.params = params
        self._model = None
        self._cp_major: int | None = None

    def _load_model(self):
        if self._model is not None:
            return self._model

        try:
            from cellpose import models  # type: ignore
        except ImportError as e:
            raise RuntimeError(
                "Cellpose is not installed. Run: pip install -r requirements.txt"
            ) from e

        self._cp_major = _cellpose_major_version()

        if self._cp_major >= 4:
            # Cellpose 4+: universal model, model_type is ignored.
            self._model = models.CellposeModel(gpu=self.params.use_gpu)
        else:
            # Cellpose 3.x: specify the model type.
            model_name = self.params.model_name
            try:
                self._model = models.CellposeModel(
                    gpu=self.params.use_gpu, model_type=model_name
                )
            except TypeError:
                # Older Cellpose API.
                self._model = models.Cellpose(
                    gpu=self.params.use_gpu, model_type=model_name
                )
            except ValueError:
                print(f"  ! Model '{model_name}' unavailable, falling back to cyto2.")
                self._model = models.CellposeModel(
                    gpu=self.params.use_gpu, model_type="cyto2"
                )

        return self._model

    def segment(
        self, image: np.ndarray, already_preprocessed: bool = False
    ) -> np.ndarray:
        """Run Cellpose on an image and return labeled masks.

        If ``already_preprocessed`` is True, the image is passed to
        Cellpose as-is (used when the caller built a composite from
        multiple channels via ``preprocess_composite``).

        Returns a 2-D labeled mask (uint16) where each cell has a unique
        integer label and background is 0.
        """
        if already_preprocessed:
            processed = image
        else:
            processed = preprocess_for_cellpose(
                image, blur_sigma=self.params.blur_sigma
            )

        model = self._load_model()

        eval_kwargs: dict = dict(
            diameter=self.params.diameter,
            flow_threshold=self.params.flow_threshold,
            cellprob_threshold=self.params.cellprob_threshold,
            min_size=self.params.min_size,
        )
        # channels was deprecated in Cellpose 4; only pass it for v3.
        if self._cp_major is not None and self._cp_major < 4:
            eval_kwargs["channels"] = [0, 0]
        if self.params.niter > 0:
            eval_kwargs["niter"] = self.params.niter

        try:
            result = model.eval(processed, **eval_kwargs)
        except TypeError:
            # Older Cellpose version doesn't recognize niter / min_size.
            eval_kwargs.pop("niter", None)
            eval_kwargs.pop("min_size", None)
            result = model.eval(processed, **eval_kwargs)

        masks = result[0]
        return np.asarray(masks, dtype=np.uint16)
