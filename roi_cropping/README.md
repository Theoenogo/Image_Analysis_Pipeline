# ROI Cropping (scaffolding — not yet implemented)

Python port of the MATLAB + ImageJ ROI cropping step. Replaces the
`H_roi_duplicate_image_all_channels_Recursive.m` + `ROI_Select_Duplicate_TIFF_Loop.ijm`
pair with a single Python script.

## Status

Scaffolding only. The Python implementation lands once the reference
MATLAB/ImageJ scripts are in [`matlab_reference/`](./matlab_reference).

## Planned behavior

- **Input**: the per-sample decon output from `preprocessing/` plus the
  user-edited ROIs from the ImageJ `Color_Merge_Automated_PreloadROIs_Adjust`
  macro:

  ```
  <Deconvoluted>/<sample>/
      gfp/gfp1_decon.tif
      cy/cy1_decon.tif
      roi_original/01.zip    ← from roi_drawing
      roi/01.zip             ← edited by the ImageJ Color_Merge macro
  ```

- **Output**: flat per-channel folders of single-cell cropped stacks, plus
  the adjusted ROI zips (re-origined to the crop-local coordinate frame):

  ```
  <Deconvoluted>/<sample>/
      Cropped/
          gfp/<pair>_<k>.tif
          cy/<pair>_<k>.tif
          roi/<pair>_<k>.zip
  ```

  The flat `Cropped/{gfp,cy}/` layout plugs straight into the next
  stage (`background_subtraction/`) and ultimately `manders_mcc/`.

- **Collapse of the old two-step flow**: the MATLAB pipeline duplicated
  each image pair N times (N = ROI count) and then ImageJ cropped each
  duplicate. In Python this intermediate duplication is unnecessary — we
  iterate over ROIs and produce one cropped cell per ROI directly.
