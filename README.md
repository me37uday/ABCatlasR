# ABCatlasR

`ABCatlasR` is an R-native interface to the
[Allen Brain Cell (ABC) Atlas](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas).
It helps you discover supported datasets, explore their metadata, select cells
that match a biological question, and retrieve expression data for downstream
single-cell analysis.

The package focuses strictly on **data access and representation**—leaving 
normalization, differential expression, integration, and trajectory analysis to 
downstream frameworks like Seurat, Monocle 3, or Bioconductor. The default output 
is a sparse **genes x cells matrix bundle** with aligned cell metadata, gene 
metadata, and provenance, with `SingleCellExperiment` and `SpatialExperiment` 
formats also supported.

## Installation

Until `ABCatlasR` is available on Bioconductor, install it from GitHub
instead:

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

## Quick start

### 1. See which datasets are available

`list_datasets()` informs about the list of datasets supported by the package, their modality,
whether spatial coordinates are available, and the expression representations
that can be retrieved.

```r
datasets <- list_datasets()
datasets
```

### 2. Load a dataset

`load_data()` downloads and assembles the **metadata only**. Expression
matrices are not downloaded at this stage, so you can first decide which cells
are relevant to your question.

Here we use the `Zhuang4` spatial mouse-brain dataset from the ABC Atlas release
validated for this version of `ABCatlasR`.

```r
abc <- load_data(
    dataset = "Zhuang4",
    release = "20261130"
)
abc
```

### 3. Explore the metadata and choose cells

Cell selection is usually the most important step. The appropriate variables
depend on the dataset and on your question of interest. For example, you might
select cells by:

- cell identity: class, subclass, supertype, cluster, or neurotransmitter;
- anatomy: brain section, region of interest, or anatomical division;
- subject information: donor, sex, age, genotype, or species;
- disease or clinical annotation, if available;
- spatial coordinates or dataset-specific quality annotations.

The full metadata table is available in `abc$cell_metadata`.

```r
names(abc$cell_metadata)
head(abc$cell_metadata)
```

For this example, we subset the metadata to **GABAergic cells** from a **single brain 
section** and take the first 100 cell labels to keep the workflow fast and 
reproducible:

```r
selected_meta <- subset(
    abc$cell_metadata,
    brain_section_label == "Zhuang-ABCA-4.002" & neurotransmitter == "GABA"
)

cells <- head(selected_meta$cell_label, 100)
length(cells)
```

Please note, the downstream functions rely on the cell IDs contained in the `cell_label` column of the metadata.

### 4. Check what will be downloaded

`plan_fetch()` is optional, but useful before working with large Atlas datasets.
It reports the source files needed for the selected cells, their cache status,
and the expected disk requirements.

```r
plan_fetch(
    abc,
    cell_ids = cells,
    data_type = "raw"
)
```

### 5. Retrieve the expression matrix

`fetch_data()` returns a matrix bundle by default. By default, **all genes available for this expression dataset are
returned**, giving a matrix suitable for downstream single-cell workflows.

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

The returned `abc_matrix` contains:

- `x$expression`: sparse `dgCMatrix`, genes x cells;
- `x$cell_metadata`: metadata aligned to matrix columns;
- `x$gene_data`: gene metadata aligned to matrix rows;
- `x$data_type`: the source expression representation;
- `x$provenance`: dataset, release, modality, and request provenance.

If your question concerns only a defined gene panel, `fetch_data()` also
accepts a `genes` argument. Gene filtering is optional; leaving `genes = NULL`
(the default) returns all available genes. Note that Atlas expression source
files are downloaded whole, so selecting fewer genes mainly reduces the size
of the returned matrix rather than the size of an uncached source download.

`ABCatlasR` returns expression values as published by the ABC Atlas and does
not normalize or transform them internally.

## Other return types

Sparse matrix is the default. When useful, the same request can return standard Bioconductor
containers by setting:

Use `return_type = "sce"` for a `SingleCellExperiment`, or
`return_type = "spatial"` for a `SpatialExperiment` in `fetch_data()`.

Raw expression is stored in a `counts` assay. Other source representations are
stored in an assay named after the data type.

## Using the matrix with Seurat or Monocle 3

`ABCatlasR` does not depend on Seurat or Monocle 3. The default raw matrix can
be passed to either framework when those packages are installed.

### Seurat

```r
seu <- SeuratObject::CreateSeuratObject(
    counts = x$expression,
    meta.data = x$cell_metadata
)
```

Seurat may sanitize some feature names internally. `x$gene_data` retains the
original ABC Atlas gene identifiers and symbols.

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

For count-based downstream workflows, retrieve `data_type = "raw"` as in the
example above.

## Downloads, caching, and reproducibility

Atlas files are downloaded on first use and cached locally. Cached files are
reused automatically, and interrupted downloads can be resumed.

Set a cache location with:

```r
options(ABCatlasR.download_base = "/path/to/abc_atlas_cache")
abc_download_dir()
```

After the required files have been cached, `offline = TRUE` prevents network
access and makes missing files an explicit error.

```r
abc <- load_data("Zhuang4", release = "20261130", offline = TRUE)

x <- fetch_data(
    abc,
    cell_ids = cells,
    data_type = "raw",
    offline = TRUE
)
```

Pinning the ABC Atlas release makes the source layout and provenance explicit.

## Learn more in the vignette

The vignette provides a more detailed guide to choosing a
dataset, selecting cells from metadata, retrieving matrices, spatial use
cases, downstream analysis, and reproducible release-aware workflows.

**[Read the ABCatlasR vignette](vignettes/ABCatlasR.Rmd)**

After installation:

```r
vignette("ABCatlasR", package = "ABCatlasR")
```


## Getting help

```r
?list_datasets
?load_data
?plan_fetch
?fetch_data
```

Please report package issues at
<https://github.com/me37uday/ABCatlasR/issues>.
