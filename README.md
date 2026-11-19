# ABCatlasR

`ABCatlasR` provides an R interface to the [Allen Brain Cell (ABC) Atlas](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas).

It lets you explore available datasets, select cells and genes, check how much data needs to be downloaded, and load expression data directly into R.

## Installation

Install the package from the `r-native` branch:

```r
install.packages("BiocManager")

BiocManager::install(
  c("rhdf5", "HDF5Array"),
  ask = FALSE,
  update = FALSE
)

install.packages("remotes")

remotes::install_github(
  "me37uday/ABCatlasR",
  ref = "r-native",
  upgrade = "never"
)
```

`rhdf5` and `HDF5Array` are Bioconductor packages used to read ABC Atlas H5AD expression files directly in R. Installing them explicitly ensures that a fresh R installation can resolve these dependencies correctly.

`SpatialExperiment` output needs the `SpatialExperiment` package, which is suggested rather than required:

```r
BiocManager::install("SpatialExperiment")
```

Then load it:

```r
library(ABCatlasR)
```

### Requirements

- **R (>= 4.4.0)** and the matching Bioconductor release.
- **No Python, compiler or system libraries.** HDF5 support comes from `rhdf5` and `Rhdf5lib`, which bundle the HDF5 library on Linux, macOS and Windows.
- **Internet access** to the public ABC Atlas bucket, `https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com`, the first time each file is used. No account or credentials are needed. Behind a proxy, set the `https_proxy` environment variable before starting R.
- **Disk space** for the local cache. Metadata files are small, but a single dataset's expression files can take several gigabytes. See [Downloads and local cache](#downloads-and-local-cache).

## Quick start

The basic workflow is: choose a dataset, load its metadata, select cells and genes while checking the download, and fetch expression.

### 1. See the available datasets

```r
list_datasets()
```

This shows the supported datasets, including which are spatial and which expression types are available.

### 2. Load metadata

```r
abc <- load_data("Zhuang4")
```

`load_data()` loads the full cell and gene metadata. It does **not** download the large expression matrices.

### 3. Select cells and genes, and check the download

```r
genes <- c("Htr7", "Sema3c", "Itga8", "Rxrg")

plan <- plan_fetch(
  abc,
  filters = list(
    brain_section_label = "Zhuang-ABCA-4.002",
    neurotransmitter = c("GABA", "Glut")
  ),
  n_cells = 500,
  seed = 1,
  genes = genes
)
```

For exact cells, use `cell_ids = my_cell_labels` instead of `filters`.

`plan_fetch()` selects the cells and genes and shows which expression files are needed, their sizes, whether they are already cached, and how much additional data must be downloaded. In an interactive session it then asks `Proceed with download? [y/n]:`. Answering `y` approves the plan; answering `n` leaves it unapproved, and the request must be planned again before fetching. `confirm = FALSE` approves the plan without asking, which is the default in non-interactive sessions such as batch jobs.

### 4. Fetch expression

```r
x <- fetch_data(plan)
```

`fetch_data()` downloads the files listed in the approved plan and reads the planned cells and genes. Downloaded files are reused from the local cache on future requests.

## Choosing an output format

### Sparse matrix

```r
x <- fetch_data(plan, return_type = "matrix")
```

The default. A list with the sparse `expression` matrix (genes x cells), `cell_metadata`, `gene_data`, `data_type` and `provenance`.

### SingleCellExperiment

```r
x <- fetch_data(plan, return_type = "sce")
```

Raw expression is stored in a `counts` assay and log2 expression in a `log2` assay.

### SpatialExperiment

```r
x <- fetch_data(plan, return_type = "spatial")
```

For spatial datasets such as `Zhuang4`, the spatial coordinates are included in the returned object. For other datasets, name the coordinate columns with `coordinate_columns`. Requires `SpatialExperiment`.

### Seurat and Monocle 3

`ABCatlasR` does not build Seurat or Monocle objects. Fetch raw expression as
a matrix bundle and convert it with the package you use:

```r
plan <- plan_fetch(abc, genes = genes, data_type = "raw")
x <- fetch_data(plan)

# Seurat (requires SeuratObject). Seurat replaces "_" in feature names.
counts <- x$expression
rownames(counts) <- gsub("_", "-", rownames(counts))
seu <- SeuratObject::CreateSeuratObject(
  counts = counts,
  meta.data = x$cell_metadata
)

# Monocle 3 (requires monocle3, installed from GitHub)
gene_metadata <- x$gene_data
gene_metadata$gene_short_name <- gene_metadata$gene_symbol
cds <- monocle3::new_cell_data_set(
  x$expression,
  cell_metadata = x$cell_metadata,
  gene_metadata = gene_metadata
)
```

Use raw expression for these count-based tools; published `log2` expression is
already transformed upstream.

## Raw and transformed expression

Use `list_datasets()` to see which representations are available.

- `"raw"` — raw counts.
- `"log2"` — log2-transformed expression published by the ABC Atlas.

If `data_type` is not specified in `plan_fetch()`, the package uses the dataset's default representation.

```r
plan <- plan_fetch(abc, genes = genes, data_type = "raw")
x <- fetch_data(plan, return_type = "matrix")
```

Some representations contain only a subset of the genes listed in the dataset metadata. `plan_fetch()` reports requested genes that are missing from the expression file.

## Downloads and local cache

Downloaded files are cached so they do not need to be downloaded again. Only the files a request needs are downloaded, and interrupted downloads resume where they stopped. Every file served from the cache is reported with a `Using cached ...; no download needed.` message. Use `plan_fetch()` to see how much a request will download before fetching.

The cache directory is chosen in this order:

1. the `download_base` argument of `load_data()`, `plan_fetch()` or `abc_cache_connect()`;
2. `options(ABCatlasR.download_base = "/path/to/abc_atlas")`;
3. the `ABC_ATLAS_CACHE` environment variable;
4. a per-user cache directory, `tools::R_user_dir("ABCatlasR", which = "cache")`.

Point the cache at a disk with enough free space before fetching expression data, for example:

```r
options(ABCatlasR.download_base = "/path/to/abc_atlas")
```

Check the active location:

```r
abc_download_dir()
```

Cached files are reused automatically. To also forbid network access, for example on a compute node without internet, use `offline = TRUE`; a file that is not cached then becomes an error instead of a download:

```r
abc <- load_data("Zhuang4", offline = TRUE)

plan <- plan_fetch(abc, genes = genes, offline = TRUE)
x <- fetch_data(plan)
```

## Complete example

```r
library(ABCatlasR)

list_datasets()

abc <- load_data("Zhuang4")

genes <- c("Htr7", "Sema3c", "Itga8", "Rxrg")

plan <- plan_fetch(
  abc,
  filters = list(
    brain_section_label = "Zhuang-ABCA-4.002",
    neurotransmitter = c("GABA", "Glut")
  ),
  n_cells = 500,
  seed = 1,
  genes = genes
)   # answer y to the prompt

x <- fetch_data(plan)

dim(x$expression)
x$provenance
```

## Getting help

```r
?list_datasets
?load_data
?plan_fetch
?fetch_data
```

For a more detailed walkthrough:

```r
vignette("ABCatlasR", package = "ABCatlasR")
```
