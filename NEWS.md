# ABCatlasR 0.99.0

## New package

* Provides an R-native, release-aware interface to Allen Brain Cell Atlas data.
* Supports metadata discovery and expression retrieval for selected cells and genes.
* Supports raw and upstream published transformed expression representations where available.
* Supports sparse matrix, SingleCellExperiment and SpatialExperiment output; the vignette and `?fetch_data` show how to convert the matrix output to Seurat and Monocle 3 objects.
* Provides local caching, resumable downloads, offline access to previously downloaded data, and an optional report of the files and disk space a request needs (`plan_fetch()`).
* `fetch_data()` checks free disk space before downloading and keeps 2 GB free on the system disk. The same check runs before every other file download, including metadata.
* A cached file whose size differs from the manifest is treated as damaged: it is downloaded again, or reported when offline.
* File paths in the manifest that point outside the cache folder are refused.
* Routes are validated against ABC Atlas release 20261130, the default. The HMBA datasets are also read from releases up to 20260711, which use earlier folder names; `load_data()` picks the folder layout that matches the requested release. When a release lacks a folder a dataset needs, `load_data()` says so and names the validated release.
* The spatial HMBA datasets offer `data_type = "raw"` (counts), `"batch_normalized"` and `"ln"` (natural log of the batch-normalized values), as published from release 20261031 on. Their cell types come from the MapMyCells results published with them.
* New datasets `HMBA_Subcortex_Marmoset` (10x Multiome) and `HMBA_spatial_Subcortex_Marmoset` (Xenium), with their own taxonomy, for 25 datasets in total.
* Preserves dataset, release, representation, and source provenance in returned objects.
