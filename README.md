# ABCatlas

`ABCatlas` provides an R interface to the [Allen Brain Cell (ABC) Atlas](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas).

It lets you explore available datasets, select cells and genes, check how much data needs to be downloaded, and load expression data directly into R.

**No Python installation is required.**

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
  "me37uday/ABCatlas",
  ref = "r-native",
  upgrade = "never"
)
```

`rhdf5` and `HDF5Array` are Bioconductor packages used to read ABC Atlas H5AD expression files directly in R. Installing them explicitly ensures that a fresh R installation can resolve these dependencies correctly. No Python installation is required.

Then load it:

```r
library(ABCatlas)
```

### Requirements

- **R (>= 4.6.0)** and the matching Bioconductor release.
- **No Python, compiler or system libraries.** HDF5 support comes from `rhdf5` and `Rhdf5lib`, which bundle the HDF5 library on Linux, macOS and Windows. `SingleCellExperiment` and `SpatialExperiment` outputs work out of the box.
- **Internet access** to the public ABC Atlas bucket, `https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com`, the first time each file is used. No account or credentials are needed. Behind a proxy, set the `https_proxy` environment variable before starting R.
- **Disk space** for the local cache. Metadata files are small, but a single dataset's expression files can take several gigabytes. See [Downloads and local cache](#downloads-and-local-cache).

## Quick start

The basic workflow is: choose a dataset, select cells, check the download, and fetch expression.

### 1. See the available datasets

```r
list_datasets()
```

This shows the supported datasets, including which are spatial and which expression types are available.

### 2. Choose cells

```r
abc <- load_data(
  "AgingMouse",
  filters = list(anatomical_division_label = "STR"),
  n_cells = 500,
  seed = 1
)
```

`load_data()` loads metadata and selects cells. It does **not** download the large expression matrices.

For exact cells:

```r
abc <- load_data(
  "AgingMouse",
  cell_ids = my_cell_labels
)
```

### 3. Check the download

```r
genes <- c("P2ry12", "Cx3cr1", "Tmem119")

plan_fetch(
  abc,
  genes = genes
)
```

`plan_fetch()` shows which expression files are needed, their sizes, whether they are already cached, and how much additional data must be downloaded. It does **not** start the download.

### 4. Fetch expression

```r
x <- fetch_data(
  abc,
  genes = genes,
  return_type = "matrix"
)
```

Required expression files are downloaded if necessary and reused from the local cache on future requests.

## Choosing an output format

### Sparse matrix

```r
x <- fetch_data(abc, genes = genes, return_type = "matrix")
```

This is the simplest option and does not require an additional analysis framework.

### SingleCellExperiment

```r
x <- fetch_data(abc, genes = genes, return_type = "sce")
```

Requires `SingleCellExperiment`.

### SpatialExperiment

```r
spatial <- load_data(
  "Zhuang4",
  n_cells = 500,
  seed = 1
)

x <- fetch_data(
  spatial,
  genes = genes,
  return_type = "spatial"
)
```

Requires `SpatialExperiment`. Spatial coordinates provided by the dataset are included in the returned object.

### Seurat and Monocle 3

`ABCatlas` does not build Seurat or Monocle objects itself, which keeps its
dependencies small. Fetch raw expression as a matrix bundle and convert it
with the package you use:

```r
x <- fetch_data(abc, genes = genes, data_type = "raw")

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

- `"raw"` — raw expression with count semantics.
- `"log2"` — expression transformed upstream by the data provider.

If `data_type` is not specified, the package uses the dataset's default representation.

```r
x <- fetch_data(
  abc,
  genes = genes,
  data_type = "raw",
  return_type = "matrix"
)
```

The package preserves the representation supplied by the ABC Atlas and does not silently convert transformed expression into counts.

Some representations may contain expression for only a subset of the genes in the broader dataset metadata. If a requested gene is not available in that expression matrix, `fetch_data()` reports it rather than silently dropping it.

## Downloads and local cache

Downloaded files are cached so they do not need to be downloaded again. Only the files a request needs are downloaded, and interrupted downloads resume where they stopped. Every file served from the cache is reported with a `Using cached ...; no download needed.` message. Use `plan_fetch()` to see how much a request will download before fetching.

The cache directory is chosen in this order:

1. the `download_base` argument of `load_data()`, `plan_fetch()`, `fetch_data()` or `abc_cache_connect()`;
2. `options(ABCatlas.download_base = "/path/to/abc_atlas")`;
3. the `ABC_ATLAS_CACHE` environment variable;
4. a per-user cache directory, `tools::R_user_dir("ABCatlas", which = "cache")`.

Point the cache at a disk with enough free space before fetching expression data, for example:

```r
options(ABCatlas.download_base = "/path/to/abc_atlas")
```

Check the active location:

```r
abc_download_dir()
```

Cached files are reused automatically. To also forbid network access, for example on a compute node without internet, use `offline = TRUE`; a file that is not cached then becomes an error instead of a download:

```r
abc <- load_data("AgingMouse", offline = TRUE)

x <- fetch_data(
  abc,
  genes = genes,
  offline = TRUE
)
```

## Complete example

```r
library(ABCatlas)

list_datasets()

abc <- load_data(
  "AgingMouse",
  filters = list(anatomical_division_label = "STR"),
  n_cells = 500,
  seed = 1
)

genes <- c("P2ry12", "Cx3cr1", "Tmem119")

plan_fetch(abc, genes = genes)

x <- fetch_data(
  abc,
  genes = genes,
  return_type = "matrix"
)

x
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
vignette("ABCatlas", package = "ABCatlas")
```
