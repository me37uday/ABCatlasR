# ABCatlasR 0.99.0

## New package

* Developed with AI assistance (Assisted-by: Anthropic Claude); all code was reviewed, tested and is maintained by the authors.
* Provides an R-native, release-aware interface to Allen Brain Cell Atlas data.
* Supports metadata discovery and expression retrieval for selected cells and genes.
* Supports raw and upstream published transformed expression representations where available.
* Supports sparse matrix, SingleCellExperiment and SpatialExperiment output; the vignette and `?fetch_data` show how to convert the matrix output to Seurat and Monocle 3 objects.
* Provides local caching, resumable downloads, offline access to previously downloaded data, and an optional report of the files and disk space a request needs (`plan_fetch()`).
* SpatialExperiment is installed with the package, so `return_type = "spatial"` works without further setup.
* `inst/conda/environment-hpc.yml` describes a conda environment with every dependency pre-built, for installing on HPC systems without compiling.
* When the ABC Atlas server cannot be reached, `load_data()`, `plan_fetch()` and `fetch_data()` continue offline with the files already downloaded, with a message, so cached data stay usable on compute nodes without internet access.
* `fetch_data()` checks free disk space before downloading and keeps 2 GB free on the system disk. The same check runs before every other file download, including metadata.
* A cached file whose size differs from the manifest is treated as damaged: it is downloaded again, or reported when offline.
* File paths in the manifest that point outside the cache folder are refused.
* Routes are validated against ABC Atlas release 20261130, the default. The HMBA datasets are also read from releases up to 20260711, which use earlier folder names; `load_data()` picks the folder layout that matches the requested release. When a release lacks a folder a dataset needs, `load_data()` says so and names the validated release.
* The spatial HMBA datasets offer `data_type = "raw"` (counts), `"batch_normalized"` and `"ln"` (natural log of the batch-normalized values), as published from release 20261031 on. Their cell types come from the MapMyCells results published with them.
* Includes the `HMBA_Subcortex_Marmoset` (10x Multiome) and `HMBA_spatial_Subcortex_Marmoset` (Xenium) datasets released in 20261130, with their own taxonomy; 25 datasets in total.
* Preserves dataset, release, representation, and source provenance in returned objects.
