# CellTrek 0.0.95

## Seurat v5 compatibility

CellTrek now runs on Seurat v5 (tested with Seurat 5.3.0 / SeuratObject 5.2.0 on R 4.3).
Verified end-to-end on VisiumV1 (brain_st_cortex + brain_sc) and VisiumV2 (stxBrain +
allen_cortex) data through traint, celltrek, scoloc (KL/DT/KD), scoexp (cc) and run_kdist.

* Replaced direct assay-slot access (`obj[[assay]]@data` / `@counts` / `@scale.data`)
  with the version-agnostic `GetAssayData()` / `SetAssayData()` accessors throughout
  `celltrek.R`, `scoexp.R` and `FindCorMarkers.R`. The old v3 idiom of moving counts into
  the data slot now uses `SetAssayData(..., layer = 'data')`.
* Made dimensional-reduction keys v5-compliant (`celltrekraw_`, `celltrek_`, `pca_`,
  `umap_`, `tsne_`).
* `CreateSeuratObject()` inputs are coerced with `as.matrix()` where needed.

## Spatial coordinates

* `traint()` now reads ST coordinates via `GetTissueCoordinates()`, supporting both
  VisiumV1 (`imagerow`/`imagecol`) and VisiumV2/FOV (`x`/`y`) image classes. Writes to
  `image@coordinates` are guarded to VisiumV1 (VisiumV2 has no such slot).
* New `traint(coord_df = )` argument: supply ST spatial coordinates directly to bypass
  differences/naming issues across image structures. Resolution priority is
  `coord_df` > auto-detect > a clear error asking for `coord_df`. Rownames must cover all
  ST cells/spots; columns named in `coord_xy` are used if present, otherwise the first two
  columns are taken as (coord_x, coord_y).

## Dependencies

* Interpolation now uses `interp::interpp` instead of `akima::interpp` (akima is not
  available on all platforms, e.g. macOS arm64). Updated DESCRIPTION and NAMESPACE.
