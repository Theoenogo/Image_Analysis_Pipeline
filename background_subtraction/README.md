# Background Subtraction (scaffolding — not yet implemented)

Python port of the ImageJ background-subtraction macro that runs between
`roi_cropping/` and `manders_mcc/`.

## Status

Scaffolding only. The Python implementation lands once the reference
ImageJ macro is in [`imagej_reference/`](./imagej_reference).

## Planned behavior

- **Input**: cropped single-cell stacks from `roi_cropping/`:

  ```
  <Deconvoluted>/<sample>/Cropped/
      gfp/<pair>_<k>.tif
      cy/<pair>_<k>.tif
      roi/<pair>_<k>.zip
  ```

- **Output**: flat per-channel folders of background-subtracted stacks,
  directly consumable by `manders_mcc/`:

  ```
  <Deconvoluted>/<sample>/Background_Subtracted/
      gfp/<pair>_<k>.tif
      cy/<pair>_<k>.tif
  ```

Exact algorithm (rolling-ball radius vs. per-cell constant subtraction
vs. something else) will match whatever the current ImageJ macro does.
