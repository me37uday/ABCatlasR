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
- **No compiler or system libraries.** HDF5 support comes from `rhdf5` and `Rhdf5lib`, which bundle the HDF5 library on Linux, macOS and Windows.
- **Internet access** to the public ABC Atlas bucket, `https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com`, the first time each file is used. No account or credentials are needed. Behind a proxy, set the `https_proxy` environment variable before starting R.
- **Disk space** for the local cache. Metadata files are small, but a single dataset's expression files can take several gigabytes. See [Downloads and local cache](#downloads-and-local-cache).

## Quick start

The basic workflow is: choose a dataset, load its metadata, choose cells and genes, optionally check the download, and fetch expression.

```r
library(ABCatlasR)

# 1. See the available datasets.
list_datasets()

# 2. Load metadata. Expression is not downloaded yet.
abc <- load_data(
  dataset = "Zhuang4",
  download_base = "~/abc_atlas_cache",
  release = "20261130"
)

# 3. Choose cells by their cell_label. Here, every cell of the dataset.
head(abc$cell_metadata)
cells <- abc$cell_metadata$cell_label

# 4. Optional: see the files these cells need, their size, cache status
#    and the disk space required. Nothing is downloaded.
plan_fetch(
  data = abc,
  cell_ids = cells,
  data_type = "raw"
)

# 5. Download what is needed and read the whole matrix into R:
#    every selected cell and, with genes = NULL, every gene.
x <- fetch_data(
  data = abc,
  cell_ids = cells,
  genes = NULL,
  data_type = "raw",
  return_type = "matrix"
)

dim(x$expression)       # genes x cells, sparse
head(x$cell_metadata)   # one row per cell
head(x$gene_data)       # one row per gene
```

The result is ready for your own downstream analysis. To fetch fewer cells or genes, pass a subset instead:

```r
# GABAergic cells only
cells <- abc$cell_metadata$cell_label[
  abc$cell_metadata$neurotransmitter %in% "GABA"
]

# A few genes, by symbol or gene identifier
x <- fetch_data(
  data = abc,
  cell_ids = cells,
  genes = c("Htr7", "Sema3c", "Itga8", "Rxrg")
)
```

### Example: a forest plot of one gene by cell group

The log2-transformed expression published by the ABC Atlas can go straight into an analysis. Here, the mean log2 expression of Htr7 in each neurotransmitter group of Zhuang4, with 95% confidence intervals, in base R:

```r
x <- fetch_data(
  data = abc,
  cell_ids = abc$cell_metadata$cell_label,
  genes = "Htr7",
  data_type = "log2"
)

expr  <- as.numeric(x$expression[1, ])  # the one requested gene
group <- x$cell_metadata$neurotransmitter
keep  <- !is.na(group) & nzchar(group)

stats <- do.call(rbind, lapply(split(expr[keep], group[keep]), function(v) {
  n <- length(v)
  half <- qt(0.975, n - 1) * sd(v) / sqrt(n)
  data.frame(n = n, mean = mean(v), lower = mean(v) - half, upper = mean(v) + half)
}))
stats <- stats[stats$n > 1, ]
stats <- stats[order(stats$mean), ]

y <- seq_len(nrow(stats))
op <- par(mar = c(4, 9, 1, 1))
plot(stats$mean, y, xlim = range(stats$lower, stats$upper), ylim = c(0.5, max(y) + 0.5),
  pch = 15, yaxt = "n", ylab = "", xlab = "Htr7, log2 expression (mean and 95% CI)")
segments(stats$lower, y, stats$upper, y)
axis(2, at = y, labels = paste0(rownames(stats), " (n = ", stats$n, ")"), las = 1)
par(op)
```

The intervals describe the spread of cells within one dataset. Cells from the same animal are not independent replicates, so compare groups across animals or datasets before drawing biological conclusions.

- `load_data()` loads the full cell and gene metadata of a dataset. The returned object records the release and the cache folder, which `plan_fetch()` and `fetch_data()` reuse.
- Without `release`, `load_data()` uses 20261130, the ABC Atlas release this version of the package is validated against. Naming the release, as above, keeps an analysis reproducible. Older releases can be requested too; release 20261031 renamed the folders of the HMBA datasets, and `load_data()` reads earlier releases with the earlier folder names.
- `cell_ids` is any vector of `cell_label` values, either all of them or a subset of the metadata. `abc$unique_values` lists the values of the main metadata columns.
- `plan_fetch()` is optional; `fetch_data()` does not need it. Without `cell_ids`, it reports the files of every cell in the dataset. Expression files are downloaded whole, so the report does not depend on the genes.
- `fetch_data()` requires `cell_ids`; to fetch every cell, use `cell_ids = abc$cell_metadata$cell_label`. With `genes = NULL`, the default, every gene in `abc$gene_data` is returned. The whole matrix of Zhuang4 is small, but for large single-cell datasets every cell and every gene can need tens of GB of memory.
- Before downloading, `fetch_data()` checks the free space of the disk holding the cache. A download that does not fit is refused, and on the disk the operating system runs from at least 2 GB must remain free afterwards, so the system keeps room for virtual memory. Files already cached need no space.

## Choosing an output format

`return_type` selects the output of `fetch_data()`.

### Sparse matrix

```r
x <- fetch_data(abc, cell_ids = cells, return_type = "matrix")
```

The default. A list with the sparse `expression` matrix (genes x cells), `cell_metadata`, `gene_data`, `data_type` and `provenance`.

### SingleCellExperiment

```r
x <- fetch_data(abc, cell_ids = cells, return_type = "sce")
```

Raw expression is stored in a `counts` assay; other data types are stored in an assay of the same name, such as `log2`.

### SpatialExperiment

```r
x <- fetch_data(abc, cell_ids = cells, return_type = "spatial")
```

Available for spatial datasets such as `Zhuang4`; the spatial coordinates are included in the returned object. `SpatialExperiment` is suggested rather than required, so install it first:

```r
BiocManager::install("SpatialExperiment")
```

### Seurat and Monocle 3

`ABCatlasR` does not build Seurat or Monocle objects. Fetch raw expression as
a matrix bundle and convert it with the package you use:

```r
x <- fetch_data(abc, cell_ids = cells, data_type = "raw")

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

Use raw expression for these count-based tools; published `log2` and `ln`
values are already transformed upstream.

## Raw and transformed expression

Use `list_datasets()` to see which representations are available.

- `"raw"` — raw counts.
- `"log2"` — log2-transformed expression published by the ABC Atlas.
- `"batch_normalized"` — counts normalized for differences between batches. Spatial HMBA datasets only.
- `"ln"` — natural log of the batch-normalized values. Spatial HMBA datasets only.

From release 20261031 on, the spatial HMBA datasets publish `"raw"`, `"batch_normalized"` and `"ln"` instead of `"log2"`.

If `data_type` is not specified, the dataset's default is used: `"raw"` where it is published, otherwise `"log2"`.

If a selected gene is not in the expression file, `fetch_data()` stops with an error that names the file.

## Downloads and local cache

Downloaded files are cached so they do not need to be downloaded again. Only the files a request needs are downloaded, and interrupted downloads resume where they stopped. Every file served from the cache is reported with a `Using cached ...; no download needed.` message. Use `plan_fetch()` to see how much a request will download before fetching.

The cache directory is chosen in this order:

1. the `download_base` argument of `load_data()` or `abc_cache_connect()`;
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

x <- fetch_data(abc, cell_ids = cells, offline = TRUE)
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
