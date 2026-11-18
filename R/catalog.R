#' List supported ABC Atlas dataset routes
#'
#' Returns the curated dataset routes supported by the native R backend.
#' This function reads the package registry only and performs no network
#' access or downloads.
#'
#' @return A data frame with one row per supported dataset route and columns
#'   describing modality, spatial status, published expression
#'   representations, and native metadata/expression routing.
#'
#' @examples
#' datasets <- list_datasets()
#' datasets[datasets$spatial, ]
#'
#' @export
list_datasets <- function() {
  routes <- .abc_datasets
  rows <- lapply(names(routes), function(dataset) {
    .abc_catalog_row(dataset, routes[[dataset]])
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("abc_dataset_catalog", "data.frame")
  out
}

# One catalog row describing a dataset route: modality, published expression
# representations and native metadata/expression directories.
.abc_catalog_row <- function(dataset, spec) {
  available <- spec$available_data_types
  if (is.null(available) || !length(available)) {
    available <- c("raw", "log2")
  }
  default <- spec$default_data_type
  if (is.null(default) || !length(default)) {
    default <- available[[1L]]
  }
  taxonomy <- if (length(spec$taxonomy)) {
    as.character(spec$taxonomy[[1L]])
  } else {
    NA_character_
  }

  data.frame(
    dataset = dataset,
    modality = as.character(spec$modality %||% "transcriptomics"),
    spatial = isTRUE(spec$spatial),
    imputed = identical(spec$modality, "imputed_spatial_transcriptomics"),
    default_data_type = as.character(default[[1L]]),
    available_data_types = paste(as.character(available), collapse = ","),
    cell_directory = .abc_spec_dir(spec$cell),
    gene_directory = .abc_spec_dir(spec$gene),
    taxonomy = taxonomy,
    expression_directories = paste(
      as.character(spec$expression_directories %||% character()),
      collapse = ","
    ),
    stringsAsFactors = FALSE
  )
}

# Directory of a registry entry such as spec$cell or spec$gene, or NA.
.abc_spec_dir <- function(entry) {
  if (!is.null(entry) && !is.null(entry$dir)) {
    as.character(entry$dir)
  } else {
    NA_character_
  }
}

#' @export
print.abc_dataset_catalog <- function(x, ...) {
  cat(
    sprintf(
      "ABC Atlas datasets: %d routes (%d spatial, %d imputed)\n\n",
      nrow(x),
      sum(x$spatial),
      sum(x$imputed)
    )
  )

  display <- data.frame(
    dataset = x$dataset,
    modality = x$modality,
    spatial = ifelse(
      x$spatial,
      "yes",
      "no"
    ),
    data = x$available_data_types,
    stringsAsFactors = FALSE
  )

  print.data.frame(
    display,
    row.names = FALSE,
    right = FALSE,
    ...
  )

  invisible(x)
}
