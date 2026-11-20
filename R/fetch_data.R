.abc_require_output_package <- function(package, return_type) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop(
            "return_type = \"", return_type, "\" needs the ", package,
            " package. Install it with BiocManager::install(\"", package,
            "\") and run the request again.",
            call. = FALSE
        )
    }
}

# Spatial output needs a dataset route declared spatial; coordinates of other
# datasets (for example UMAP embeddings) are not anatomical locations.
.abc_require_spatial <- function(cell_metadata) {
    if (!isTRUE(attr(cell_metadata, "abc_spatial"))) {
        dataset <- attr(cell_metadata, "abc_dataset") %||% "This dataset"
        stop(
            dataset, " has no spatial coordinates, so return_type = ",
            "\"spatial\" is not available. Use \"matrix\" or \"sce\"; ",
            "list_datasets() shows which datasets are spatial.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.abc_spatial_coordinates <- function(metadata) {
    chosen <- .abc_default_spatial_columns(metadata)

    coordinates <- as.matrix(metadata[, chosen, drop = FALSE])
    storage.mode(coordinates) <- "double"
    if (anyNA(coordinates) || any(!is.finite(coordinates))) {
        stop(
            "Spatial coordinates must be finite and non-missing.",
            call. = FALSE
        )
    }

    rownames(coordinates) <- as.character(metadata$cell_label)
    colnames(coordinates) <- chosen
    coordinates
}

# First published coordinate pair present in the metadata, in order of
# preference: slab (mm), rotated (um), experiment, then plain x/y.
.abc_default_spatial_columns <- function(metadata) {
    candidates <- list(
        c("x_slab_mm", "y_slab_mm"),
        c("x_rotated_um", "y_rotated_um"),
        c("x_experiment", "y_experiment"),
        c("x", "y")
    )
    present <- vapply(
        candidates, function(x) all(x %in% names(metadata)), logical(1)
    )
    if (!any(present)) {
        stop(
            "No supported spatial coordinate pair was found in the cell ",
            "metadata.",
            call. = FALSE
        )
    }

    candidates[[which(present)[[1L]]]]
}

.abc_output_provenance <- function(
        metadata, gene_data, data_type,
        release = NULL) {
    list(
        dataset = attr(metadata, "abc_dataset"),
        release = release,
        data_type = data_type,
        modality = attr(metadata, "abc_modality"),
        spatial = isTRUE(attr(metadata, "abc_spatial")),
        expression_directories = attr(
            metadata,
            "abc_expression_directories"
        ),
        n_cells = nrow(metadata),
        n_genes = nrow(gene_data)
    )
}

.abc_as_dgc_matrix <- function(expression) {
    if (inherits(expression, "dgCMatrix")) {
        return(expression)
    }

    expression <- methods::as(
        expression,
        "generalMatrix"
    )

    methods::as(
        expression,
        "CsparseMatrix"
    ) |>
        methods::as("dgCMatrix")
}

.abc_materialization_gene_ids <- function(gene_data) {
    gene_ids <- as.character(
        gene_data$gene_identifier
    )

    if (anyNA(gene_ids) ||
        any(!nzchar(gene_ids)) ||
        anyDuplicated(gene_ids)) {
        stop(
            "gene_data$gene_identifier must be a non-missing unique key.",
            call. = FALSE
        )
    }

    gene_ids
}

.abc_validate_materialization <- function(
        expression,
        cell_metadata,
        gene_data
) {
    if (ncol(expression) != nrow(cell_metadata)) {
        stop(
            "Expression columns do not match the selected cell metadata.",
            call. = FALSE
        )
    }

    if (nrow(expression) != nrow(gene_data)) {
        stop(
            "Expression rows do not match the selected gene metadata.",
            call. = FALSE
        )
    }

    if (!identical(
        as.character(colnames(expression)),
        as.character(cell_metadata$cell_label)
    )) {
        stop(
            "Expression columns are not aligned to cell_metadata$cell_label.",
            call. = FALSE
        )
    }

    gene_ids <- .abc_materialization_gene_ids(
        gene_data
    )

    if (!identical(
        as.character(rownames(expression)),
        gene_ids
    )) {
        stop(
            "Expression rows are not aligned to gene_data$gene_identifier.",
            call. = FALSE
        )
    }

    gene_ids
}

.abc_default_materialization_provenance <- function(
        cell_metadata,
        gene_data,
        data_type
) {
    list(
        data_type = data_type,
        modality = attr(
            cell_metadata,
            "abc_modality"
        ),
        spatial = isTRUE(
            attr(cell_metadata, "abc_spatial")
        ),
        n_cells = nrow(cell_metadata),
        n_genes = nrow(gene_data)
    )
}

.abc_prepare_materialization <- function(
        expression,
        cell_metadata,
        gene_data,
        data_type,
        provenance = NULL
) {
    expression <- .abc_as_dgc_matrix(
        expression
    )

    gene_ids <- .abc_validate_materialization(
        expression,
        cell_metadata,
        gene_data
    )

    rownames(cell_metadata) <- as.character(
        cell_metadata$cell_label
    )
    rownames(gene_data) <- gene_ids

    if (is.null(provenance)) {
        provenance <- .abc_default_materialization_provenance(
            cell_metadata,
            gene_data,
            data_type
        )
    }

    list(
        expression = expression,
        cell_metadata = cell_metadata,
        gene_data = gene_data,
        provenance = provenance
    )
}

.abc_materialize_matrix <- function(
        expression,
        cell_metadata,
        gene_data,
        data_type,
        provenance
) {
    structure(
        list(
            expression = expression,
            cell_metadata = cell_metadata,
            gene_data = gene_data,
            data_type = data_type,
            provenance = provenance
        ),
        class = c("abc_matrix", "list")
    )
}

.abc_bioc_assays <- function(expression, data_type) {
    assay_name <- if (data_type == "raw") {
        "counts"
    } else {
        "log2"
    }

    stats::setNames(
        list(expression),
        assay_name
    )
}

.abc_materialize_sce <- function(
        expression,
        cell_metadata,
        gene_data,
        data_type,
        provenance
) {
    .abc_require_output_package(
        "SingleCellExperiment",
        "sce"
    )
    .abc_require_output_package(
        "S4Vectors",
        "sce"
    )

    SingleCellExperiment::SingleCellExperiment(
        assays = .abc_bioc_assays(
            expression,
            data_type
        ),
        rowData = S4Vectors::DataFrame(gene_data),
        colData = S4Vectors::DataFrame(cell_metadata),
        metadata = list(
            ABCatlasR = provenance
        )
    )
}

.abc_materialize_spatial <- function(
        expression,
        cell_metadata,
        gene_data,
        data_type,
        provenance
) {
    .abc_require_output_package(
        "SpatialExperiment",
        "spatial"
    )
    .abc_require_spatial(cell_metadata)

    coordinates <- .abc_spatial_coordinates(cell_metadata)

    provenance$coordinate_columns <- colnames(
        coordinates
    )

    SpatialExperiment::SpatialExperiment(
        assays = .abc_bioc_assays(
            expression,
            data_type
        ),
        rowData = S4Vectors::DataFrame(gene_data),
        colData = S4Vectors::DataFrame(cell_metadata),
        spatialCoords = coordinates,
        metadata = list(
            ABCatlasR = provenance
        )
    )
}

.abc_dispatch_materialization <- function(
        prepared,
        data_type,
        return_type
) {
    switch(return_type,
        matrix = .abc_materialize_matrix(
            prepared$expression,
            prepared$cell_metadata,
            prepared$gene_data,
            data_type,
            prepared$provenance
        ),
        sce = .abc_materialize_sce(
            prepared$expression,
            prepared$cell_metadata,
            prepared$gene_data,
            data_type,
            prepared$provenance
        ),
        spatial = .abc_materialize_spatial(
            prepared$expression,
            prepared$cell_metadata,
            prepared$gene_data,
            data_type,
            prepared$provenance
        )
    )
}

.abc_materialize <- function(
        expression,
        cell_metadata,
        gene_data,
        data_type,
        return_type = "matrix",
        provenance = NULL
) {
    return_type <- match.arg(
        return_type,
        c("matrix", "sce", "spatial")
    )

    prepared <- .abc_prepare_materialization(
        expression = expression,
        cell_metadata = cell_metadata,
        gene_data = gene_data,
        data_type = data_type,
        provenance = provenance
    )

    .abc_dispatch_materialization(
        prepared = prepared,
        data_type = data_type,
        return_type = return_type
    )
}


# Resolve the cache, release, representation and cells of a request on `data`.
# plan_fetch() and fetch_data() each call this; neither uses the other.
.abc_request_context <- function(data, cell_ids, data_type, offline) {
    if (!inherits(data, "abc_data")) {
        stop(
            "data must be an abc_data object returned by load_data().",
            call. = FALSE
        )
    }

    download_base <- abc_download_dir(data$download_base)
    release <- data$release %||% .abc_default_release
    abc_cache <- abc_cache_connect(
        download_base,
        offline = offline, release = release
    )
    cells <- .abc_select_cells(data$cell_metadata, cell_ids)

    list(
        cells = cells,
        download_base = download_base,
        release = release,
        data_type = .abc_resolve_data_type(data_type, data$cell_metadata),
        abc_cache = abc_cache
    )
}

# Expression files needed for `cells`, with the bytes still to download.
.abc_request_files <- function(context) {
    files <- .abc_expression_files(
        context$abc_cache, context$cells,
        data_type = context$data_type
    )
    files$additional_download <- ifelse(
        files$cached, 0, pmax(files$size - files$partial_bytes, 0)
    )
    files
}

.abc_fetch_expression <- function(context, files, gene_data, gene_sel) {
    total_download <- sum(files$additional_download)
    storage <- .abc_storage_check(context$download_base, total_download)

    if (!storage$sufficient) {
        stop(.abc_storage_message(storage), call. = FALSE)
    }
    if (total_download > 0 && is.na(storage$available)) {
        message(
            "Free disk space could not be determined; downloading ",
            .abc_format_size(total_download), "."
        )
    }

    .abc_expression_matrix(
        abc_cache = context$abc_cache, cells = context$cells,
        gene_data = gene_data, gene_sel = gene_sel,
        data_type = context$data_type
    )
}

# Seurat and Monocle objects are built by the user from the matrix bundle
# (see "Converting to other formats" in ?fetch_data). Asking for one of them
# points there instead of giving match.arg()'s generic error.
.abc_match_return_type <- function(return_type) {
    if (length(return_type) == 1L &&
        tolower(return_type) %in% c("monocle", "monocle3", "seurat")) {
        stop(
            "return_type = \"", return_type, "\" is not built in. ",
            "Use return_type = \"matrix\" (or \"sce\") and convert the ",
            "result; see \"Converting to other formats\" in ?fetch_data.",
            call. = FALSE
        )
    }
    match.arg(return_type, c("matrix", "sce", "spatial"))
}

# Validate the requested output before any expression file is resolved or
# downloaded, so a request that cannot be returned fails without a transfer.
# The materializers repeat these checks as a second line of defence.
.abc_validate_output_request <- function(return_type, cell_metadata) {
    if (identical(return_type, "sce")) {
        .abc_require_output_package("SingleCellExperiment", "sce")
    }

    if (identical(return_type, "spatial")) {
        .abc_require_output_package("SpatialExperiment", "spatial")
        .abc_require_spatial(cell_metadata)
        .abc_spatial_coordinates(cell_metadata)
    }

    invisible(TRUE)
}

#' Fetch ABC Atlas expression for selected cells and genes
#'
#' Downloads the expression files the selected cells need, if they are not
#' already cached, and reads the selected cells and genes into the requested
#' R representation. [plan_fetch()] reports the same files and sizes
#' beforehand, but is not needed to fetch.
#'
#' Before downloading, the free space of the disk holding the cache is
#' checked. A download that does not fit is refused. On the disk the
#' operating system runs from, at least 2 GB must also remain free after
#' the download, so the system keeps room for virtual memory. No space is
#' checked when every file is already cached. The same check runs before
#' any other file is downloaded, including the metadata read by
#' [load_data()].
#'
#' @param data An \code{abc_data} object returned by [load_data()]. It
#'   provides the metadata, release and cache directory.
#' @param cell_ids Character vector of cell labels (the
#'   \code{cell_label} column of \code{data$cell_metadata}) to fetch. The
#'   columns of the result follow this order. To fetch every cell, use
#'   \code{data$cell_metadata$cell_label}.
#' @param genes Character vector of genes to include. Each entry is matched
#'   first against \code{gene_identifier}, then against the published source
#'   identifier, then against \code{gene_symbol}; a symbol shared by several
#'   genes is reported as ambiguous. If \code{NULL} (default), all genes of
#'   the expression file are returned.
#' @param data_type Published expression representation, \code{"raw"} or
#'   \code{"log2"}. If \code{NULL} (default), the dataset's default is used:
#'   \code{"raw"} where it is published, otherwise \code{"log2"}. Raw
#'   expression keeps count semantics; published log2 expression is never
#'   presented as counts.
#' @param return_type Output format: \code{"matrix"} (default) for the sparse
#'   expression bundle, \code{"sce"} for a SingleCellExperiment, or
#'   \code{"spatial"} for a SpatialExperiment with the dataset's spatial
#'   coordinates. \code{"spatial"} needs the SpatialExperiment package and a
#'   spatial dataset (see [list_datasets()]). For Seurat or Monocle 3
#'   objects, see "Converting to other formats" below.
#' @param offline If \code{TRUE}, use only locally cached files and never
#'   contact S3; a file that is not cached is an error instead of a download.
#' @return The requested expression representation: a sparse matrix bundle
#'   (a list with \code{expression}, \code{cell_metadata},
#'   \code{gene_data}, \code{data_type} and \code{provenance}), a
#'   SingleCellExperiment, or a SpatialExperiment.
#' @section Converting to other formats:
#' Seurat and Monocle 3 objects are not built by this package. Fetch raw
#' expression as a matrix bundle and convert it with the package you use.
#' Both tools expect count data, so fetch with \code{data_type = "raw"}.
#'
#' ```r
#' x <- fetch_data(abc, cell_ids = cells, genes = genes, data_type = "raw")
#'
#' # Seurat (package SeuratObject). Seurat replaces "_" in feature names.
#' counts <- x$expression
#' rownames(counts) <- gsub("_", "-", rownames(counts))
#' seu <- SeuratObject::CreateSeuratObject(
#'     counts = counts,
#'     meta.data = x$cell_metadata
#' )
#'
#' # Monocle 3 (package monocle3, installed from GitHub).
#' gene_metadata <- x$gene_data
#' gene_metadata$gene_short_name <- gene_metadata$gene_symbol
#' cds <- monocle3::new_cell_data_set(
#'     x$expression,
#'     cell_metadata = x$cell_metadata,
#'     gene_metadata = gene_metadata
#' )
#' ```
#' @examples
#' # A small synthetic dataset shipped with the package.
#' example_dir <- system.file("extdata", "ex", package = "ABCatlasR")
#' abc <- readRDS(file.path(example_dir, "example-data.rds"))
#' abc$download_base <- example_dir
#'
#' cells <- c("c4", "c1")
#' x <- fetch_data(
#'     abc,
#'     cell_ids = cells,
#'     genes = c("G1", "G4"),
#'     offline = TRUE
#' )
#' x$expression
#' @export
fetch_data <- function(
        data,
        cell_ids,
        genes = NULL,
        data_type = NULL,
        return_type = c("matrix", "sce", "spatial"),
        offline = FALSE) {
    if (missing(cell_ids)) {
        stop(
            "cell_ids is required. To fetch every cell, use ",
            "cell_ids = data$cell_metadata$cell_label.",
            call. = FALSE
        )
    }
    if (is.null(cell_ids)) {
        stop("cell_ids must contain at least one cell label.", call. = FALSE)
    }
    return_type <- .abc_match_return_type(return_type)

    context <- .abc_request_context(data, cell_ids, data_type, offline)
    selected_genes <- .abc_resolve_genes(data$gene_data, genes)
    .abc_validate_output_request(return_type, context$cells)

    expression <- .abc_fetch_expression(
        context,
        files = .abc_request_files(context),
        gene_data = data$gene_data,
        gene_sel = selected_genes$.abc_gene_index
    )
    selected_genes$.abc_gene_index <- NULL
    provenance <- .abc_output_provenance(
        context$cells, selected_genes,
        data_type = context$data_type, release = context$release
    )

    .abc_materialize(
        expression = expression, cell_metadata = context$cells,
        gene_data = selected_genes, data_type = context$data_type,
        return_type = return_type, provenance = provenance
    )
}

#' Report the expression files a request needs
#'
#' Lists the expression files that [fetch_data()] would need for the selected
#' cells, with their sizes, whether they are already cached, and how much
#' would still be downloaded, and checks that against the free disk space.
#' Nothing is downloaded. Calling \code{plan_fetch()} is optional:
#' [fetch_data()] does not use its result.
#'
#' Expression files are downloaded whole, so the sizes do not depend on the
#' genes that will be fetched. They depend on the cells, because the cells of
#' some datasets are spread over several files, and on the representation.
#'
#' @param data An \code{abc_data} object returned by [load_data()].
#' @param cell_ids Optional character vector of cell labels. If \code{NULL}
#'   (default), every cell of the dataset is included.
#' @param data_type Published expression representation, \code{"raw"} or
#'   \code{"log2"}. If \code{NULL} (default), the dataset's default is used:
#'   \code{"raw"} where it is published, otherwise \code{"log2"}.
#' @param offline If \code{TRUE}, use only locally cached files and never
#'   contact S3.
#' @return An \code{abc_fetch_plan} object, printed as a report of the
#'   expression files, their sizes and cache status, and the storage check.
#' @examples
#' # A small synthetic dataset shipped with the package.
#' example_dir <- system.file("extdata", "ex", package = "ABCatlasR")
#' abc <- readRDS(file.path(example_dir, "example-data.rds"))
#' abc$download_base <- example_dir
#'
#' plan_fetch(abc, offline = TRUE)
#'
#' # Exact cells: only the files holding them are listed.
#' plan_fetch(abc, cell_ids = c("c1", "c2"), offline = TRUE)
#' @export
plan_fetch <- function(
        data, cell_ids = NULL, data_type = NULL, offline = FALSE) {
    context <- .abc_request_context(data, cell_ids, data_type, offline)
    files <- .abc_request_files(context)
    total_download <- sum(files$additional_download)

    structure(
        list(
            dataset = attr(context$cells, "abc_dataset"),
            release = context$abc_cache$version,
            data_type = context$data_type,
            n_cells = nrow(context$cells),
            files = files,
            total_size = sum(files$size),
            total_partial = sum(files$partial_bytes),
            total_download = total_download,
            storage = .abc_storage_check(context$download_base, total_download)
        ),
        class = "abc_fetch_plan"
    )
}

#' @export
print.abc_fetch_plan <- function(x, ...) {
    cat(
        "<abc_fetch_plan>\n",
        "  dataset:        ", x$dataset %||% "", "\n",
        "  release:        ", x$release, "\n",
        "  representation: ", x$data_type, "\n",
        "  cells:          ", format(x$n_cells, big.mark = ","), "\n",
        "\nExpression files:\n",
        sep = ""
    )
    for (i in seq_len(nrow(x$files))) {
        .abc_print_plan_file(x$files[i, ])
    }

    storage <- x$storage
    space <- if (is.na(storage$available)) {
        "unable to determine"
    } else {
        .abc_format_size(storage$available)
    }
    check <- if (x$total_download == 0) {
        "PASS (nothing to download)"
    } else if (is.na(storage$available)) {
        "unable to determine"
    } else if (storage$sufficient) {
        "PASS"
    } else {
        "INSUFFICIENT"
    }

    cat(
        "\nTotals:\n",
        "  selected-file size:   ", .abc_format_size(x$total_size), "\n",
        "  resumable partial:    ", .abc_format_size(x$total_partial), "\n",
        "  additional download:  ", .abc_format_size(x$total_download), "\n",
        "  available disk space: ", space, "\n",
        "  required free space:  ", .abc_format_size(storage$required), "\n",
        "  storage check:        ", check, "\n",
        "\nNo expression files were downloaded.\n",
        sep = ""
    )
    invisible(x)
}

# One expression-file entry of a printed fetch plan.
.abc_print_plan_file <- function(f) {
    cat(
        "  ", f$dataset_label, "/", f$feature_matrix_label, "/", f$data_type,
        " [", f$n_cells, " cells]\n",
        "    size:                ", .abc_format_size(f$size), "\n",
        "    cached:              ", if (f$cached) "yes" else "no", "\n",
        "    resumable partial:   ", .abc_format_size(f$partial_bytes), "\n",
        "    additional download: ", .abc_format_size(f$additional_download),
        "\n",
        sep = ""
    )
}

# Resolve an expression representation from dataset metadata.
.abc_resolve_data_type <- function(data_type, metadata) {
    default <- attr(metadata, "abc_default_data_type")
    available <- attr(metadata, "abc_available_data_types")

    if (is.null(default) || !length(default)) {
        default <- "raw"
    }
    if (is.null(available) || !length(available)) {
        available <- c("raw", "log2")
    }

    default <- as.character(default[[1L]])
    available <- as.character(available)

    if (is.null(data_type) || !length(data_type) ||
        !nzchar(as.character(data_type[[1L]]))) {
        data_type <- default
    }

    if (length(data_type) != 1L || is.na(data_type)) {
        stop("data_type must be one of: 'raw', 'log2'.", call. = FALSE)
    }

    data_type <- as.character(data_type)

    if (!data_type %in% c("raw", "log2")) {
        stop("data_type must be one of: 'raw', 'log2'.", call. = FALSE)
    }

    if (!data_type %in% available) {
        dataset <- attr(metadata, "abc_dataset")
        if (is.null(dataset) || !length(dataset)) {
            dataset <- "This dataset"
        }

        stop(
            dataset, " does not publish data_type='", data_type,
            "'. Available data_type values: ",
            paste(sprintf("'%s'", available), collapse = ", "),
            "; default: ", default, ".",
            call. = FALSE
        )
    }

    data_type
}

# Resolve requested genes against stable identifiers first, then source
# identifiers when available, and finally symbols. Each request must identify
# exactly one gene; absent genes are never silently replaced by zero rows.
.abc_resolve_gene_request <- function(
        query,
        ids,
        source_ids,
        symbols
) {
    matches <- which(ids == query)

    if (!length(matches) && !is.null(source_ids)) {
        matches <- which(
            !is.na(source_ids) &
                source_ids == query
        )
    }

    if (!length(matches) && !is.null(symbols)) {
        matches <- which(
            !is.na(symbols) &
                symbols == query
        )
    }

    if (!length(matches)) {
        stop(
            "Gene '",
            query,
            "' was not found.",
            call. = FALSE
        )
    }

    if (length(matches) > 1L) {
        choices <- utils::head(
            ids[matches],
            8L
        )

        stop(
            "Gene '",
            query,
            "' is ambiguous. Please use one of these gene identifiers: ",
            paste(choices, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    matches[[1L]]
}

.abc_validate_gene_resolution_inputs <- function(gene_data, genes) {
    if (!"gene_identifier" %in% names(gene_data)) {
        stop("'gene_identifier' column not found in gene_data.", call. = FALSE)
    }

    ids <- as.character(gene_data$gene_identifier)
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
        stop(
            "'gene_identifier' must be a non-missing unique key in gene_data.",
            call. = FALSE
        )
    }

    query <- if (is.null(genes)) NULL else as.character(genes)
    if (!is.null(query) && (!length(query) || anyNA(query) ||
        any(!nzchar(query)) || anyDuplicated(query))) {
        stop(
            "genes must be a nonempty vector without duplicate requests.",
            call. = FALSE
        )
    }

    list(ids = ids, query = query)
}

.abc_finalize_gene_resolution <- function(
        gene_data,
        selected
) {
    out <- gene_data[
        selected, ,
        drop = FALSE
    ]

    if (
        anyDuplicated(
            as.character(out$gene_identifier)
        )
    ) {
        stop(
            "Multiple requests resolve to the same gene identifier.",
            call. = FALSE
        )
    }

    out$.abc_gene_index <- selected
    rownames(out) <- NULL
    out
}

.abc_resolve_genes <- function(gene_data, genes = NULL) {
    inputs <- .abc_validate_gene_resolution_inputs(gene_data, genes)

    if (is.null(inputs$query)) {
        out <- gene_data
        out$.abc_gene_index <- seq_len(nrow(out))
        rownames(out) <- NULL
        return(out)
    }

    selected <- vapply(
        inputs$query, .abc_resolve_gene_request, integer(1L),
        ids = inputs$ids,
        source_ids = .abc_column_or_null(gene_data, "source_gene_identifier"),
        symbols = .abc_column_or_null(gene_data, "gene_symbol"),
        USE.NAMES = FALSE
    )

    .abc_finalize_gene_resolution(gene_data = gene_data, selected = selected)
}

# A column of `data` as character, or NULL when the column is absent.
.abc_column_or_null <- function(data, column) {
    if (column %in% names(data)) as.character(data[[column]]) else NULL
}

# Select cells consistently for metadata loading, planning and fetching.
.abc_validate_cell_id_selection <- function(
        metadata,
        cell_ids
) {
    if (!is.character(cell_ids)) {
        stop(
            "cell_ids must be a character vector.",
            call. = FALSE
        )
    }

    if (
        !length(cell_ids) ||
            anyNA(cell_ids) ||
            any(!nzchar(cell_ids)) ||
            anyDuplicated(cell_ids)
    ) {
        stop(
            "cell_ids must be non-missing, nonempty and contain no duplicates.",
            call. = FALSE
        )
    }

    if (!"cell_label" %in% names(metadata)) {
        stop(
            "Cell metadata is missing required column 'cell_label'.",
            call. = FALSE
        )
    }

    labels <- as.character(
        metadata$cell_label
    )

    if (
        anyNA(labels) ||
            anyDuplicated(labels)
    ) {
        stop(
            "cell_label must be a non-missing unique key before ",
            "cell_ids selection.",
            call. = FALSE
        )
    }

    labels
}

.abc_match_cell_ids <- function(
        metadata,
        cell_ids,
        labels
) {
    pos <- match(
        cell_ids,
        labels
    )

    if (anyNA(pos)) {
        missing <- cell_ids[is.na(pos)]

        stop(
            length(missing),
            " requested cell",
            if (length(missing) == 1L) "" else "s",
            " not found in this dataset: ",
            paste(
                utils::head(missing, 3L),
                collapse = ", "
            ),
            if (length(missing) > 3L) ", ..." else "",
            call. = FALSE
        )
    }

    metadata[
        pos, ,
        drop = FALSE
    ]
}

# Select the cells named in `cell_ids`, in that order; NULL selects every
# cell of the dataset.
.abc_select_cells <- function(metadata, cell_ids = NULL) {
    out <- if (is.null(cell_ids)) {
        metadata
    } else {
        .abc_match_cell_ids(
            metadata = metadata,
            cell_ids = cell_ids,
            labels = .abc_validate_cell_id_selection(metadata, cell_ids)
        )
    }

    if (!nrow(out)) {
        stop("The dataset has no cells.", call. = FALSE)
    }

    rownames(out) <- as.character(out$cell_label)

    .abc_status_done(
        sprintf(
            "Selected %s cells.",
            format(nrow(out), big.mark = ",")
        )
    )

    out
}

# Resolve the expression files needed for `cells` without downloading them.
.abc_validate_expression_matrix_labels <- function(cells) {
    if (!"feature_matrix_label" %in% names(cells)) {
        stop(
            "Cell metadata is missing required column: feature_matrix_label.",
            call. = FALSE
        )
    }

    matrix_labels <- as.character(
        cells$feature_matrix_label
    )

    if (
        anyNA(matrix_labels) ||
            any(!nzchar(matrix_labels))
    ) {
        stop(
            "Selected cells contain a missing feature_matrix_label.",
            call. = FALSE
        )
    }

    matrix_labels
}

.abc_apply_expression_matrix_map <- function(
        matrix_labels,
        cells
) {
    matrix_map <- attr(
        cells,
        "abc_expression_matrix_map"
    )

    if (is.null(matrix_map) || !length(matrix_map)) {
        return(matrix_labels)
    }

    map_names <- names(matrix_map)

    matrix_map <- stats::setNames(
        as.character(matrix_map),
        map_names
    )

    if (
        is.null(map_names) ||
            anyNA(map_names) ||
            any(!nzchar(map_names)) ||
            anyDuplicated(map_names) ||
            anyNA(matrix_map) ||
            any(!nzchar(matrix_map))
    ) {
        stop(
            "Invalid ABC expression-matrix mapping provenance.",
            call. = FALSE
        )
    }

    mapped <- unname(
        matrix_map[matrix_labels]
    )

    use_map <- !is.na(mapped)
    matrix_labels[use_map] <- mapped[use_map]

    matrix_labels
}

.abc_expression_matrix_labels <- function(cells) {
    matrix_labels <- .abc_validate_expression_matrix_labels(
        cells
    )

    .abc_apply_expression_matrix_map(
        matrix_labels = matrix_labels,
        cells = cells
    )
}

.abc_legacy_expression_route_pairs <- function(
        cells,
        matrix_labels
) {
    if (!"dataset_label" %in% names(cells)) {
        stop(
            "Cell metadata has no ABC expression-route provenance ",
            "and is missing ",
            "dataset_label. Use metadata returned by load_data(), ",
            "or provide a ",
            "dataset_label column for legacy routing.",
            call. = FALSE
        )
    }

    dataset_labels <- as.character(
        cells$dataset_label
    )

    if (
        anyNA(dataset_labels) ||
            any(!nzchar(dataset_labels))
    ) {
        stop(
            "Selected cells contain a missing dataset_label.",
            call. = FALSE
        )
    }

    unique(
        data.frame(
            directory = dataset_labels,
            feature_matrix_label = matrix_labels,
            stringsAsFactors = FALSE
        )
    )
}

.abc_expression_route_inventory <- function(
        abc_cache,
        expression_dirs,
        data_type
) {
    expression_dirs <- unique(
        as.character(expression_dirs)
    )

    inventory <- do.call(
        rbind,
        lapply(
            expression_dirs,
            .abc_expression_inventory,
            cache = abc_cache
        )
    )

    inventory[
        inventory$data_type == data_type, ,
        drop = FALSE
    ]
}

.abc_match_expression_route_label <- function(
        inventory,
        label,
        expression_dirs,
        data_type
) {
    matches <- inventory[
        inventory$feature_matrix_label == label,
        c("directory", "feature_matrix_label"),
        drop = FALSE
    ]

    if (nrow(matches) == 0L) {
        stop(
            "Cannot route feature_matrix_label='",
            label,
            "' for data_type='",
            data_type,
            "' within this dataset's expression directories: ",
            paste(expression_dirs, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    if (nrow(matches) > 1L) {
        choices <- paste(
            paste0(
                matches$directory,
                "/",
                matches$feature_matrix_label
            ),
            collapse = ", "
        )

        stop(
            "Ambiguous expression routing for feature_matrix_label='",
            label,
            "' and data_type='",
            data_type,
            "': ",
            choices,
            ".",
            call. = FALSE
        )
    }

    matches
}

.abc_inventory_expression_route_pairs <- function(
        abc_cache,
        expression_dirs,
        matrix_labels,
        data_type
) {
    expression_dirs <- unique(
        as.character(expression_dirs)
    )

    inventory <- .abc_expression_route_inventory(
        abc_cache = abc_cache,
        expression_dirs = expression_dirs,
        data_type = data_type
    )

    labels <- unique(matrix_labels)

    pairs <- lapply(
        labels,
        .abc_match_expression_route_label,
        inventory = inventory,
        expression_dirs = expression_dirs,
        data_type = data_type
    )

    do.call(rbind, pairs)
}

.abc_expression_route_pairs <- function(
        abc_cache,
        cells,
        matrix_labels,
        data_type
) {
    expression_dirs <- attr(
        cells,
        "abc_expression_directories"
    )

    if (is.null(expression_dirs) || !length(expression_dirs)) {
        return(
            .abc_legacy_expression_route_pairs(
                cells = cells,
                matrix_labels = matrix_labels
            )
        )
    }

    .abc_inventory_expression_route_pairs(
        abc_cache = abc_cache,
        expression_dirs = expression_dirs,
        matrix_labels = matrix_labels,
        data_type = data_type
    )
}

.abc_expression_files <- function(abc_cache, cells, data_type = "raw") {
    matrix_labels <- .abc_expression_matrix_labels(cells)
    pairs <- .abc_expression_route_pairs(
        abc_cache = abc_cache, cells = cells,
        matrix_labels = matrix_labels, data_type = data_type
    )
    pairs <- pairs[
        order(pairs$directory, pairs$feature_matrix_label, method = "radix"), ,
        drop = FALSE
    ]

    out <- vector("list", nrow(pairs))
    for (k in seq_len(nrow(pairs))) {
        out[[k]] <- .abc_expression_file_row(
            abc_cache,
            directory = as.character(pairs$directory[k]),
            matrix_file = as.character(pairs$feature_matrix_label[k]),
            data_type = data_type,
            matrix_labels = matrix_labels
        )
    }
    do.call(rbind, out)
}

# Size and local cache state (complete file and resumable .part bytes) of one
# expression matrix file, as one row of the fetch plan's file table.
.abc_expression_file_row <- function(
        abc_cache, directory, matrix_file,
        data_type, matrix_labels) {
    attrs <- .abc_file_attributes(
        abc_cache, directory, paste0(matrix_file, "/", data_type)
    )
    part <- paste0(attrs$local_path, ".part")
    partial <- if (file.exists(part)) as.numeric(file.info(part)$size) else 0

    data.frame(
        dataset_label = directory,
        feature_matrix_label = matrix_file,
        data_type = data_type,
        n_cells = sum(matrix_labels == matrix_file),
        size = as.numeric(attrs$size),
        relative_path = attrs$relative_path,
        local_path = attrs$local_path,
        cached = file.exists(attrs$local_path) && .abc_size_matches(attrs),
        partial_bytes = partial,
        stringsAsFactors = FALSE
    )
}

# Open the central H5AD expression matrix without realizing it in memory.
#
# ABC Atlas H5AD files currently publish /X either as an AnnData sparse-matrix
# group (for example csr_matrix) or as a dense HDF5 dataset with
# encoding-type="array". Return a genes x cells delayed matrix in both cases.
.abc_h5ad_expression_layout <- function(path) {
    info <- rhdf5::h5ls(
        path,
        recursive = FALSE
    )

    entry <- info[
        info$group == "/" &
            info$name == "X", ,
        drop = FALSE
    ]

    if (nrow(entry) != 1L) {
        stop(
            basename(path),
            ": expected exactly one /X object in H5AD file.",
            call. = FALSE
        )
    }

    attrs <- rhdf5::h5readAttributes(
        path,
        "/X"
    )

    encoding <- attrs[["encoding-type"]]

    if (is.null(encoding) || !length(encoding)) {
        encoding <- NA_character_
    } else {
        encoding <- as.character(
            encoding[[1L]]
        )
    }

    list(
        otype = as.character(entry$otype[[1L]]),
        encoding = encoding
    )
}

.abc_h5ad_expression_from_layout <- function(
        path,
        layout
) {
    otype <- layout$otype
    encoding <- layout$encoding

    if (
        identical(otype, "H5I_GROUP") &&
            encoding %in% c(
                "csr_matrix",
                "csc_matrix"
            )
    ) {
        return(
            HDF5Array::H5SparseMatrix(
                path,
                "/X"
            )
        )
    }

    if (
        identical(otype, "H5I_DATASET") &&
            (
                is.na(encoding) ||
                    encoding %in% c(
                        "array",
                        "dense"
                    )
            )
    ) {
        return(
            HDF5Array::HDF5Array(
                path,
                "/X"
            )
        )
    }

    stop(
        basename(path),
        ": unsupported /X H5AD layout (otype=",
        otype,
        ", encoding-type=",
        if (is.na(encoding)) "<missing>" else encoding,
        ").",
        call. = FALSE
    )
}

.abc_open_h5ad_expression <- function(path) {
    layout <- .abc_h5ad_expression_layout(
        path
    )

    .abc_h5ad_expression_from_layout(
        path = path,
        layout = layout
    )
}

.abc_h5ad_expression <- function(path, n_genes = NULL, n_cells = NULL) {
    X <- .abc_open_h5ad_expression(
        path
    )

    if (
        !is.null(n_genes) &&
            !is.null(n_cells)
    ) {
        expected <- c(
            as.integer(n_genes),
            as.integer(n_cells)
        )

        observed <- as.integer(dim(X))

        if (!identical(observed, expected)) {
            stop(
                basename(path),
                ": /X dimensions are ",
                paste(observed, collapse = " x "),
                " but var/obs imply ",
                paste(expected, collapse = " x "),
                ". Automatic transposition is intentionally disabled.",
                call. = FALSE
            )
        }
    }

    X
}

# Read the requested genes x cells blocks from one H5AD expression file and
# return them as triplet pieces (row, column, value). Cells absent from the
# file contribute no entries.
.abc_expression_file_blocks <- function(
        path, cells, gene_data, gene_sel,
        block_size) {
    var_index <- .h5ad_index(path, "var")
    obs_index <- .h5ad_index(path, "obs")
    gene_pos <- .abc_gene_positions(
        gene_data, gene_sel, var_index, basename(path)
    )

    file_pos <- match(rownames(cells), obs_index)
    wanted <- which(!is.na(file_pos))
    wanted <- wanted[order(file_pos[wanted])]

    X <- .abc_h5ad_expression(
        path,
        n_genes = length(var_index), n_cells = length(obs_index)
    )
    blocks <- split(wanted, ceiling(seq_along(wanted) / block_size))
    .abc_read_expression_blocks(X, gene_pos, file_pos, blocks)
}

# Realize each block of cells for the selected genes and keep its non-zero
# entries. `blocks` holds positions in `cells`; `file_pos` maps them to columns
# of the file.
.abc_read_expression_blocks <- function(X, gene_pos, file_pos, blocks) {
    rows <- cols <- vals <- list()

    for (block in blocks) {
        m <- methods::as(
            X[gene_pos, file_pos[block], drop = FALSE],
            "dgCMatrix"
        )
        nz <- Matrix::summary(m)
        rows[[length(rows) + 1L]] <- nz$i
        cols[[length(cols) + 1L]] <- block[nz$j]
        vals[[length(vals) + 1L]] <- nz$x
    }

    list(rows = rows, cols = cols, vals = vals)
}

# Sparse genes x cells expression matrix (dgCMatrix) for `cells` (rownames are
# cell labels) and the rows of `gene_data` selected by `gene_sel`.
#
# Cells are grouped by the expression file that holds them; each file is
# downloaded if needed and opened as an on-disk delayed matrix. Both sparse
# H5AD groups and dense /X datasets are supported, and only requested
# gene/cell blocks are realized. Rows follow the requested gene order and are
# named by gene_identifier; columns follow the order of `cells`. Cells not
# found in any file are all-zero columns.
.abc_expression_matrix <- function(
        abc_cache, cells, gene_data, gene_sel,
        data_type = "raw", block_size = 50000L) {
    if (!requireNamespace("HDF5Array", quietly = TRUE)) {
        stop(
            "Package 'HDF5Array' is required to read expression matrices.",
            call. = FALSE
        )
    }

    gene_sel <- as.integer(gene_sel)
    gene_filtered <- gene_data[gene_sel, , drop = FALSE]
    cell_ids <- rownames(cells)
    files <- .abc_expression_files(abc_cache, cells, data_type = data_type)

    rows <- cols <- vals <- list()
    for (k in seq_len(nrow(files))) {
        pieces <- .abc_read_expression_file(
            abc_cache, files[k, ], cells, gene_data, gene_sel, data_type,
            block_size
        )
        rows <- c(rows, pieces$rows)
        cols <- c(cols, pieces$cols)
        vals <- c(vals, pieces$vals)
    }

    out <- Matrix::sparseMatrix(
        i = unlist(rows, use.names = FALSE) %||% integer(),
        j = unlist(cols, use.names = FALSE) %||% integer(),
        x = unlist(vals, use.names = FALSE) %||% numeric(),
        dims = c(length(gene_sel), length(cell_ids)),
        dimnames = list(as.character(gene_filtered$gene_identifier), cell_ids)
    )

    .abc_status_done(sprintf(
        "Expression ready: %s genes x %s cells.",
        format(length(gene_sel), big.mark = ","),
        format(length(cell_ids), big.mark = ",")
    ))
    out
}

# Resolve one expression file (downloading it if needed), then read its
# requested blocks under a single-line spinner. Download/cache resolution
# happens before the spinner starts so that curl's real download progress
# display is never overwritten.
.abc_read_expression_file <- function(
        abc_cache, entry, cells, gene_data,
        gene_sel, data_type, block_size) {
    matrix_file <- entry$feature_matrix_label
    path <- abc_cache$get_file_path(
        directory = entry$dataset_label,
        file_name = paste0(matrix_file, "/", data_type)
    )

    spinner <- .abc_spinner_start(
        paste0("Loading expression data: ", matrix_file)
    )
    completed <- FALSE
    tryCatch(
        {
            pieces <- .abc_expression_file_blocks(
                path = path, cells = cells, gene_data = gene_data,
                gene_sel = gene_sel, block_size = block_size
            )
            completed <- TRUE
        },
        finally = .abc_spinner_stop(
            spinner,
            if (completed) paste0("Loaded expression data: ", matrix_file, ".")
        )
    )

    pieces
}

.abc_h5ad_source_gene_identifier <- function(x) {
    sub(
        "^(NCBIGene:[^:]+):(?:NC|NW)_[^:]+$",
        "\\1",
        as.character(x)
    )
}

# Positions in an expression file's var index for selected gene metadata rows.
# Canonical gene_identifier is preferred. For normalized duplicate source IDs,
# occurrence metadata permits one-to-one alignment only when the H5AD contains
# exactly the expected number of source-ID occurrences. No aggregation or
# arbitrary duplicate selection is allowed.
#
# Some ABC Atlas files store malformed identifiers for a few genes: in
# Zeng-Aging-Mouse-10Xv3, genes with a duplicated symbol have a printed pandas
# Series embedding the real identifier. The existing positional fallback for
# that published-data defect is retained after ordinary identifier matching.
.abc_source_gene_position <- function(gene_data, row, var, var_source) {
    source <- as.character(gene_data$source_gene_identifier[[row]])
    occurrence <- as.integer(gene_data$source_gene_occurrence[[row]])
    group_size <- as.integer(gene_data$source_gene_group_size[[row]])

    if (is.na(source) || !nzchar(source) ||
        is.na(occurrence) || occurrence < 1L ||
        is.na(group_size) || group_size < 1L) {
        return(NA_integer_)
    }

    hits <- which(var == source)
    if (!length(hits)) {
        hits <- which(var_source == source)
    }

    if (length(hits) == 1L && group_size == 1L) {
        return(hits[[1L]])
    }
    if (group_size > 1L && length(hits) == group_size &&
        occurrence <= length(hits)) {
        return(hits[[occurrence]])
    }

    NA_integer_
}

.abc_resolve_source_gene_positions <- function(
        gene_data,
        sel,
        var,
        pos
) {
    required <- c(
        "source_gene_identifier",
        "source_gene_occurrence",
        "source_gene_group_size"
    )

    if (!all(required %in% names(gene_data))) {
        return(pos)
    }

    unresolved <- which(is.na(pos))

    if (!length(unresolved)) {
        return(pos)
    }

    var_source <- .abc_h5ad_source_gene_identifier(
        var
    )

    for (j in unresolved) {
        row <- sel[[j]]

        resolved <- .abc_source_gene_position(
            gene_data = gene_data,
            row = row,
            var = var,
            var_source = var_source
        )

        if (!is.na(resolved)) {
            pos[[j]] <- resolved
        }
    }

    pos
}

.abc_gene_positions <- function(gene_data, sel, var, file) {
    sel <- as.integer(sel)
    ids <- as.character(gene_data$gene_identifier %||% rownames(gene_data))

    pos <- match(ids[sel], var)
    if (!anyNA(pos)) {
        return(pos)
    }

    pos <- .abc_resolve_source_gene_positions(
        gene_data = gene_data, sel = sel, var = var, pos = pos
    )
    if (!anyNA(pos)) {
        return(pos)
    }

    .abc_positional_gene_fallback(
        ids, sel, var, file,
        n_missing = sum(is.na(pos))
    )
}

# Preserve the existing AgingMouse compatibility fallback: when every
# mismatching var entry embeds the expected identifier, align by position.
# Otherwise the selected genes are genuinely missing from the file.
.abc_positional_gene_fallback <- function(ids, sel, var, file, n_missing) {
    if (length(var) == length(ids)) {
        bad <- which(var != ids)

        if (length(bad) &&
            all(mapply(grepl, ids[bad], var[bad],
                MoreArgs = list(fixed = TRUE)
            ))) {
            warning(
                length(bad), " gene identifiers in ", file,
                " are malformed; genes were aligned by position.",
                call. = FALSE
            )
            return(sel)
        }
    }

    stop(
        n_missing, " of the selected genes are missing from ", file,
        "; is gene_data from the same dataset?",
        call. = FALSE
    )
}

.h5ad_index <- function(path, group) {
    index_col <- rhdf5::h5readAttributes(path, paste0("/", group))[["_index"]]
    if (is.null(index_col)) {
        stop(
            "Cannot find the ", group, " index in ", basename(path),
            " (unsupported h5ad layout)."
        )
    }
    as.character(rhdf5::h5read(path, paste0("/", group, "/", index_col)))
}
