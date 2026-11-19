# ABCatlasR 0.99.0

## New package

* Provides an R-native, release-aware interface to Allen Brain Cell Atlas data.
* Supports metadata discovery, cell selection, and expression-data retrieval.
* Supports raw and upstream published transformed expression representations where available.
* Supports sparse matrix, SingleCellExperiment and SpatialExperiment output; the vignette and `?fetch_data` show how to convert the matrix output to Seurat and Monocle 3 objects.
* Provides local caching, resumable downloads, offline access to previously downloaded data, and download planning.
* Preserves dataset, release, representation, and source provenance in returned objects.
