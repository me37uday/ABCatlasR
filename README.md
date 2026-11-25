# ABCatlasR

`ABCatlasR` is an R-native interface to the
[Allen Brain Cell (ABC) Atlas](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas).
It provides a simple workflow for discovering supported datasets, exploring
cell and gene metadata, selecting cells that match a biological question, and
retrieving expression data for downstream single-cell analysis.

The default output is a sparse **genes x cells** matrix together with aligned
cell metadata, gene metadata, and provenance. The same request can also be
returned as a `SingleCellExperiment` or `SpatialExperiment`.

## Installation

Upon availability in Bioconductor, install with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("ABCatlasR")
```

Until then, install the development version from GitHub:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

options(repos = BiocManager::repositories())

remotes::install_github(
    "me37uday/ABCatlasR"
)
```

Load the package:

```r
library(ABCatlasR)
```

For installation on HPC systems with Conda/Bioconda, see the
[package vignette](vignettes/ABCatlasR.Rmd).

## Quick start

### 1. Choose where Atlas files are stored

ABC Atlas files are downloaded when they are first needed and then reused from
local storage. Set the download directory before starting the analysis.

```r
options(ABCatlasR.download_base = "/path/to/abc_atlas_cache")
abc_download_dir()
```

### 2. See which datasets are available

`list_datasets()` lists the datasets supported by the installed version of
`ABCatlasR`, together with their modality, spatial status, and available
expression representations.

```r
datasets <- list_datasets()
datasets
```

### 3. Load a dataset

`load_data()` loads the cell and gene metadata for a dataset. Expression data
are not retrieved at this stage, so the metadata can first be used to define
the cells relevant to the biological question. The examples below use Zhuang4,
a spatial mouse-brain resource.

```r
abc <- load_data(
    dataset = "Zhuang4",
    release = "20261130"
)
abc
```

### 4. Explore the metadata and select cells

The available annotations depend on the dataset. Start by inspecting the cell
metadata and the curated values provided in `unique_values`.

```r
names(abc$cell_metadata)
head(abc$cell_metadata)
```

For example, select **GABAergic cells** from **one brain section**:

```r
selected_meta <- subset(
    abc$cell_metadata,
    brain_section_label == "Zhuang-ABCA-4.002" &
        neurotransmitter == "GABA"
)

cells <- selected_meta$cell_label
length(cells)
```

`cell_label` contains the cell identifiers used by `plan_fetch()` and
`fetch_data()`.

### 5. Inspect the request

`plan_fetch()` shows which expression files are required, whether they are
already available locally, and the storage associated with the request.
It is optional, but recommended before larger requests.

```r
plan_fetch(
    abc,
    cell_ids = cells,
    data_type = "raw"
)
```

### 6. Retrieve the expression matrix

`fetch_data()` returns an `abc_matrix` by default. When `genes` is not
specified, all genes available for the selected expression representation are
returned.

```r
x <- fetch_data(
    abc,
    cell_ids = cells,
    data_type = "raw"
)
```

Inspect the result:

```r
x
dim(x$expression)
class(x$expression)
head(x$cell_metadata)
head(x$gene_data)
```

The main components are:

- `x$expression`: sparse `dgCMatrix`, genes x cells;
- `x$cell_metadata`: metadata aligned to matrix columns;
- `x$gene_data`: gene metadata aligned to matrix rows;
- `x$data_type`: expression representation returned by the Atlas;
- `x$provenance`: dataset, release, modality, and request information.

`ABCatlasR` returns the expression representation published by the ABC Atlas;
it does not normalize or transform the values internally.

## Other return types

The sparse matrix bundle is the default output. The same request can be
returned in standard Bioconductor containers:

```r
sce <- fetch_data(
    abc,
    cell_ids = cells,
    data_type = "raw",
    return_type = "sce"
)
```

For datasets with spatial coordinates:

```r
spe <- fetch_data(
    abc,
    cell_ids = cells,
    data_type = "raw",
    return_type = "spatial"
)
```

Use `list_datasets()` to identify spatial datasets and the expression
representations available for each route.

## Using the matrix with Seurat or Monocle 3

`ABCatlasR` does not depend on Seurat or Monocle 3, but the default matrix can
be passed directly to either framework when those packages are installed.

### Seurat

```r
seu <- SeuratObject::CreateSeuratObject(
    counts = x$expression,
    meta.data = x$cell_metadata
)
```

### Monocle 3

```r
gene_metadata <- x$gene_data
gene_metadata$gene_short_name <- gene_metadata$gene_symbol

cds <- monocle3::new_cell_data_set(
    x$expression,
    cell_metadata = x$cell_metadata,
    gene_metadata = gene_metadata
)
```

For count-based downstream workflows, retrieve `data_type = "raw"` when raw
counts are available for the selected dataset.

## Learn more

The vignette provides a more complete description of the datasets and
functions in `ABCatlasR`, including dataset selection, metadata exploration,
expression representations, return types, release-aware caching, offline use,
and HPC installation.

**[Read the ABCatlasR vignette](vignettes/ABCatlasR.Rmd)**

After installation:

```r
vignette("ABCatlasR", package = "ABCatlasR")
```
