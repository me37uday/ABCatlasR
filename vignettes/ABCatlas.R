## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
vignette_source_root <- normalizePath(file.path(dirname(knitr::current_input()), ".."), mustWork = FALSE)

example_dir <- system.file(
  "extdata", "ex",
  package = "ABCatlas"
)
if (!nzchar(example_dir)) example_dir <- file.path(vignette_source_root, "inst", "extdata", "ex")
example_data <- readRDS(
  file.path(example_dir, "example-data.rds")
)

vignette_snapshot_path <- system.file(
  "extdata",
  "vignette",
  "AgingMouse-STR-1000-plan.rds",
  package = "ABCatlas"
 )
if (!nzchar(vignette_snapshot_path)) vignette_snapshot_path <- file.path(vignette_source_root, "inst", "extdata", "vignette", "AgingMouse-STR-1000-plan.rds")
vignette_snapshot <- readRDS(vignette_snapshot_path)

vignette_plan_cache <- system.file(
  "extdata",
  "vignette",
  "AgingMouse-plan-cache",
  package = "ABCatlas"
)
if (!nzchar(vignette_plan_cache)) vignette_plan_cache <- file.path(vignette_source_root, "inst", "extdata", "vignette", "AgingMouse-plan-cache")

## ----library------------------------------------------------------------------
library(ABCatlas)

## ----list-datasets------------------------------------------------------------
datasets <- list_datasets()
datasets

## ----spatial-datasets---------------------------------------------------------
datasets[datasets$spatial, ]

## ----load-data, eval=FALSE----------------------------------------------------
# abc <- load_data(
#   "AgingMouse"
# )
#
# abc
# abc$cell_metadata
# abc$gene_data
# abc$unique_values

## ----cell-selection-----------------------------------------------------------
abc_filtered <- load_data(
  "AgingMouse",
  filters = list(
    anatomical_division_label = "STR"
  )
)


## ----cell-id-selection, eval=FALSE--------------------------------------------
# abc_exact <- load_data(
#   "AgingMouse",
#   cell_ids = my_cell_labels
# )

## ----cell-sampling------------------------------------------------------------
abc_sampled <- load_data(
  "AgingMouse",
  filters = list(
    anatomical_division_label = "STR"
  ),
  n_cells = 1000,
  seed = 42
)


## ----plan-fetch---------------------------------------------------------------
genes <- c(
  "P2ry12",
  "Cx3cr1",
  "Tmem119",
  "Apoe",
  "Cst3",
  "Mrc1",
  "Trem2",
  "Hexb",
  "Tyrobp",
  "Sall1"
)

plan <- plan_fetch(
  metadata = vignette_snapshot$cells,
  gene_data = vignette_snapshot$genes,
  genes = vignette_snapshot$genes_requested,
  download_base = vignette_plan_cache,
  offline = TRUE,
  release = vignette_snapshot$release,
  confirm = FALSE
)

plan

## ----inspect-matrix-result, eval=FALSE----------------------------------------
# x <- fetch_data(
#   plan = plan,
#   return_type = "matrix"
# )
# x$expression
# x$cell_metadata
# x$gene_data
# x$provenance

## ----fetch-monocle, eval=FALSE------------------------------------------------
# x <- fetch_data(
#   plan = plan,
#   return_type = "monocle"
# )

## ----fetch-sce, eval=FALSE----------------------------------------------------
# x <- fetch_data(
#   plan = plan,
#   return_type = "sce"
# )

## ----spatial-output, eval=FALSE-----------------------------------------------
# spatial_abc <- load_data(
#   "MERFISH",
#   n_cells = 1000,
#   seed = 1
# )
#
# x <- fetch_data(
#   spatial_abc,
#   genes = genes,
#   return_type = "spatial"
# )
#
# SpatialExperiment::spatialCoords(x)
#
# # For routes without canonical spatial coordinates, verified columns
# # can be supplied explicitly.
# y <- fetch_data(
#   abc,
#   genes = genes,
#   return_type = "spatial",
#   coordinate_columns = c(
#     "x_experiment",
#     "y_experiment"
#   )
# )

## ----data-types---------------------------------------------------------------
as.data.frame(list_datasets())[
  ,
  c(
    "dataset",
    "default_data_type",
    "available_data_types",
    "imputed"
  )
]

## ----gene-identifiers, eval=FALSE---------------------------------------------
# head(
#   abc$gene_data[
#     ,
#     c(
#       "gene_identifier",
#       "gene_symbol"
#     )
#   ]
# )
#
# fetch_data(
#   abc,
#   genes = c(
#     "P2ry12",
#     "Cx3cr1"
#   ),
#   return_type = "matrix"
# )

## ----offline-workflow, eval=FALSE---------------------------------------------
# abc <- load_data(
#   "AgingMouse",
#   offline = TRUE
# )
#
# x <- fetch_data(
#   abc,
#   genes = genes,
#   return_type = "matrix",
#   offline = TRUE
# )

## ----session-info-------------------------------------------------------------
sessionInfo()
