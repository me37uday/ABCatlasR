# Resolve either an abc_data object or the legacy metadata/gene_data pair.
.abc_fetch_inputs <- function(
        data = NULL, metadata = NULL, gene_data = NULL,
        download_base = NULL, release = NULL) {
    if (!is.null(data)) {
        if (!inherits(data, "abc_data")) {
            stop(
                "data must be an abc_data object returned by load_data().",
                call. = FALSE
            )
        }

        if (!is.null(metadata) || !is.null(gene_data)) {
            stop(
                "Supply either data or metadata/gene_data, not both.",
                call. = FALSE
            )
        }

        metadata <- data$cell_metadata
        gene_data <- data$gene_data

        if (is.null(download_base)) {
            download_base <- data$download_base
        }

        if (is.null(release)) {
            release <- data$release
        }
    }

    if (is.null(metadata) || is.null(gene_data)) {
        stop(
            "Supply an abc_data object via data, ",
            "or both metadata and gene_data.",
            call. = FALSE
        )
    }

    if (!is.data.frame(metadata)) {
        stop("metadata must be a data.frame.", call. = FALSE)
    }

    if (!is.data.frame(gene_data)) {
        stop("gene_data must be a data.frame.", call. = FALSE)
    }

    list(
        metadata = metadata,
        gene_data = gene_data,
        download_base = download_base,
        release = release
    )
}


.abc_require_output_package <- function(packages, return_type) {
    for (package in packages) {
        if (!requireNamespace(package, quietly = TRUE)) {
            stop(
                "return_type='", return_type,
                "' requires package '", package,
                "'. Install it to use this output format.",
                call. = FALSE
            )
        }
    }
}

.abc_spatial_coordinates <- function(metadata, coordinate_columns = NULL) {
    chosen <- if (is.null(coordinate_columns)) {
        .abc_default_spatial_columns(metadata)
    } else {
        .abc_check_spatial_columns(metadata, coordinate_columns)
    }

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

# Explicit coordinate columns must be two names present in the metadata.
.abc_check_spatial_columns <- function(metadata, coordinate_columns) {
    if (!is.character(coordinate_columns) ||
        length(coordinate_columns) != 2L ||
        anyNA(coordinate_columns) ||
        any(!nzchar(coordinate_columns))) {
        stop(
            "coordinate_columns must contain exactly two non-missing ",
            "column names.",
            call. = FALSE
        )
    }

    missing <- setdiff(coordinate_columns, names(metadata))
    if (length(missing)) {
        stop(
            "Spatial coordinate column",
            if (length(missing) == 1L) " " else "s ",
            paste(missing, collapse = ", "), " not found in cell metadata.",
            call. = FALSE
        )
    }

    coordinate_columns
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
            "No supported spatial coordinate pair was found. ",
            "Provide coordinate_columns explicitly.",
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

.abc_materialize_monocle <- function(
        expression,
        cell_metadata,
        gene_data,
        data_type,
        provenance
) {
    if (!identical(data_type, "raw")) {
        stop(
            "return_type = \"monocle\" requires raw expression data; ",
            "the published log2 representation cannot be used as ",
            "a count-based Monocle 2 CellDataSet.",
            call. = FALSE
        )
    }

    .abc_require_output_package(c("monocle", "Biobase"), "monocle")

    monocle_gene_data <- gene_data

    if (!"gene_short_name" %in% names(monocle_gene_data)) {
        if ("gene_symbol" %in% names(monocle_gene_data)) {
            monocle_gene_data$gene_short_name <- as.character(
                monocle_gene_data$gene_symbol
            )
        } else {
            monocle_gene_data$gene_short_name <- rownames(
                monocle_gene_data
            )
        }
    }

    cds <- monocle::newCellDataSet(
        expression,
        phenoData = Biobase::AnnotatedDataFrame(
            data = cell_metadata
        ),
        featureData = Biobase::AnnotatedDataFrame(
            data = monocle_gene_data
        )
    )

    experiment_data <- Biobase::experimentData(cds)

    experiment_data@other$ABCatlas <- provenance

    Biobase::experimentData(cds) <- experiment_data

    cds
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
            ABCatlas = provenance
        )
    )
}

.abc_materialize_spatial <- function(
        expression,
        cell_metadata,
        gene_data,
        data_type,
        coordinate_columns,
        provenance,
        spatial
) {
    .abc_require_output_package(
        "SpatialExperiment",
        "spatial"
    )
    .abc_require_output_package(
        "S4Vectors",
        "spatial"
    )

    if (is.null(coordinate_columns) && !isTRUE(spatial)) {
        stop(
            "No declared spatial modality. Supply verified coordinate_columns ",
            "explicitly; dimensionality-reduction coordinates are not ",
            "anatomical ",
            "spatial locations.",
            call. = FALSE
        )
    }

    coordinates <- .abc_spatial_coordinates(
        cell_metadata,
        coordinate_columns = coordinate_columns
    )

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
            ABCatlas = provenance
        )
    )
}

.abc_dispatch_materialization <- function(
        prepared,
        data_type,
        return_type,
        coordinate_columns,
        spatial
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
            coordinate_columns,
            prepared$provenance,
            spatial
        ),
        monocle = .abc_materialize_monocle(
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
        coordinate_columns = NULL,
        provenance = NULL,
        spatial = FALSE
) {
    return_type <- match.arg(
        return_type,
        c("matrix", "sce", "spatial", "monocle")
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
        return_type = return_type,
        coordinate_columns = coordinate_columns,
        spatial = spatial
    )
}


.abc_fetch_context_from_plan <- function(plan) {
    if (!inherits(plan, "abc_fetch_plan")) {
        stop(
            "plan must be an abc_fetch_plan object.",
            call. = FALSE
        )
    }

    if (!isTRUE(plan$approved)) {
        stop(
            "The fetch plan was not approved. ",
            "Run plan_fetch() and approve the download first.",
            call. = FALSE
        )
    }

    abc_cache <- abc_cache_connect(
        plan$.download_base,
        offline = plan$.offline,
        release = plan$release
    )

    list(
        filtered_meta = plan$.metadata,
        gene_data = plan$.gene_data,
        selected_gene_data = plan$.selected_gene_data,
        download_base = plan$.download_base,
        release = plan$release,
        data_type = plan$data_type,
        abc_cache = abc_cache
    )
}

.abc_fetch_context_from_inputs <- function(
        data, metadata, filters,
        gene_data, genes, data_type,
        offline, release, download_base) {
    # Backward compatibility for the historical positional download_base.
    if (!is.null(data) && !inherits(data, "abc_data") &&
        is.character(data) && length(data) == 1L && is.null(download_base)) {
        download_base <- data
        data <- NULL
    }

    inputs <- .abc_fetch_inputs(
        data = data, metadata = metadata, gene_data = gene_data,
        download_base = download_base, release = release
    )
    download_base <- abc_download_dir(inputs$download_base)
    release <- inputs$release %||% .abc_default_release
    abc_cache <- abc_cache_connect(
        download_base,
        offline = offline, release = release
    )
    prepared <- .abc_prepare_fetch(
        inputs$metadata, filters, inputs$gene_data, genes
    )

    list(
        filtered_meta = prepared$metadata,
        gene_data = inputs$gene_data,
        selected_gene_data = prepared$gene_data,
        download_base = download_base,
        release = release,
        data_type = .abc_resolve_data_type(data_type, prepared$metadata),
        abc_cache = abc_cache
    )
}

.abc_fetch_expression <- function(
        abc_cache, filtered_meta, gene_data,
        selected_gene_data, download_base,
        data_type) {
    files <- .abc_expression_files(abc_cache, filtered_meta,
        data_type = data_type
    )
    additional_download <- ifelse(
        files$cached, 0, pmax(files$size - files$partial_bytes, 0)
    )
    total_download <- sum(additional_download)
    available_space <- .abc_available_space(download_base)

    if (total_download > 0 && !is.na(available_space) &&
        available_space < total_download) {
        stop(
            "Insufficient disk space for this download. ",
            "Required: ", .abc_format_size(total_download),
            "; available: ", .abc_format_size(available_space),
            "; shortfall: ", .abc_format_size(total_download - available_space),
            ". No expression files were downloaded.",
            call. = FALSE
        )
    }

    .abc_expression_matrix(
        abc_cache = abc_cache, cells = filtered_meta, gene_data = gene_data,
        gene_sel = selected_gene_data$.abc_gene_index, data_type = data_type
    )
}

# Validate the requested output representation before any expression file is
# resolved or downloaded, so incompatible requests fail without a transfer.
# The materializers repeat these checks as a second line of defence.
.abc_validate_output_request <- function(
        return_type,
        data_type,
        cell_metadata,
        coordinate_columns
) {
    if (identical(return_type, "monocle") && !identical(data_type, "raw")) {
        stop(
            "return_type = \"monocle\" requires raw expression data; ",
            "the published log2 representation cannot be used as ",
            "a count-based Monocle 2 CellDataSet.",
            call. = FALSE
        )
    }

    if (identical(return_type, "monocle")) {
        .abc_require_output_package(c("monocle", "Biobase"), "monocle")
    }

    if (identical(return_type, "spatial")) {
        if (is.null(coordinate_columns) &&
            !isTRUE(attr(cell_metadata, "abc_spatial"))) {
            stop(
                "No declared spatial modality. ",
                "Supply verified coordinate_columns ",
                "explicitly; dimensionality-reduction coordinates are not ",
                "anatomical ",
                "spatial locations.",
                call. = FALSE
            )
        }

        .abc_spatial_coordinates(
            cell_metadata,
            coordinate_columns = coordinate_columns
        )
    }

    invisible(TRUE)
}

# A fetch plan already fixes the cells, genes, release, cache directory and
# representation. Arguments that would otherwise be silently ignored are
# rejected so the user knows the plan, not those arguments, defines the request.
.abc_reject_plan_overrides <- function(
        data,
        metadata,
        filters,
        gene_data,
        genes,
        data_type,
        offline,
        release,
        download_base
) {
    supplied <- c(
        data = !is.null(data),
        metadata = !is.null(metadata),
        filters = length(filters) > 0L,
        gene_data = !is.null(gene_data),
        genes = !is.null(genes),
        data_type = !is.null(data_type),
        offline = !identical(offline, FALSE),
        release = !is.null(release),
        download_base = !is.null(download_base)
    )

    if (any(supplied)) {
        stop(
            "When plan is supplied, the plan defines the request. Remove: ",
            paste(names(supplied)[supplied], collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    invisible(TRUE)
}

#' Fetch selected ABC Atlas expression into an R analysis representation
#'
#' @param data Optional \code{abc_data} object returned by \code{load_data()}.
#'   This is the preferred interface. When supplied, its metadata, gene data,
#'   cache directory, and release are reused automatically.
#' @param download_base Directory where ABC Atlas data is downloaded and cached.
#'   If `NULL` (default), it is resolved by [abc_download_dir()] (option,
#'   `ABC_ATLAS_CACHE` env var, or a per-user cache directory).
#' @param metadata A data.frame containing ABC atlas cell metadata (e.g., the
#'   `cell_metadata` element returned by [load_data()]).
#' @param filters A named list of column = value pairs used to filter the
#'   metadata.
#' @param gene_data A data.frame with 'gene_identifier' and 'gene_symbol'
#'   columns, typically the `gene_data` element returned by [load_data()].
#' @param genes Character vector of genes to include. Each entry is matched
#'   first against \code{gene_identifier}, then against the published source
#'   identifier, then against \code{gene_symbol}; a symbol shared by several
#'   genes is reported as ambiguous. If \code{NULL}, all genes are used.
#' @param return_type Output format: \code{"matrix"} for the primary sparse
#'   expression bundle, \code{"sce"} for a SingleCellExperiment,
#'   \code{"spatial"} for a SpatialExperiment, or \code{"monocle"} for
#'   a Monocle 2 CellDataSet. The default is \code{"matrix"}.
#' @param coordinate_columns Optional length-two character vector naming the
#'   spatial x/y columns. Used only for \code{return_type = "spatial"}.
#' @param data_type Published expression representation (`"raw"` or `"log2"`).
#'   If `NULL`, the dataset route's declared default is used. Raw expression
#'   retains count semantics. Upstream log2 expression remains transformed
#'   data and is never presented as counts.
#' @param offline If `TRUE`, use only locally cached expression data and never
#'   contact S3 (see [abc_cache_connect()]).
#' @param release Optional ABC Atlas release version in \code{YYYYMMDD} form.
#'   Defaults to the release validated against this package. Advanced users may
#'   supply another release explicitly.
#' @param plan Optional \code{abc_fetch_plan} object returned by
#'   [plan_fetch()] and approved interactively. When supplied, the cells,
#'   genes, release, cache directory and representation stored in the plan
#'   are used, and the download is not confirmed again.
#' @return The requested expression representation: a sparse matrix bundle,
#'   SingleCellExperiment, SpatialExperiment, or Monocle 2 CellDataSet.
#' @examples
#' example_dir <- system.file(
#'     "extdata", "ex",
#'     package = "ABCatlas"
#' )
#' example_data <- readRDS(
#'     file.path(example_dir, "example-data.rds")
#' )
#'
#' x <- fetch_data(
#'     download_base = example_dir,
#'     metadata = example_data$cells,
#'     gene_data = example_data$genes,
#'     genes = c("G1", "G4"),
#'     return_type = "matrix",
#'     offline = TRUE,
#'     release = "20250101"
#' )
#' x$expression
#' @export
fetch_data <- function(
        data = NULL, metadata = NULL, filters = list(),
        gene_data = NULL, genes = NULL,
        return_type = c("matrix", "sce", "spatial", "monocle"),
        coordinate_columns = NULL, data_type = NULL,
        offline = FALSE, release = NULL, download_base = NULL,
        plan = NULL) {
    return_type <- match.arg(return_type)
    context <- .abc_fetch_context(
        plan = plan, data = data, metadata = metadata, filters = filters,
        gene_data = gene_data, genes = genes, data_type = data_type,
        offline = offline, release = release, download_base = download_base
    )
    .abc_validate_output_request(
        return_type = return_type, data_type = context$data_type,
        cell_metadata = context$filtered_meta,
        coordinate_columns = coordinate_columns
    )

    expression <- .abc_fetch_expression(
        abc_cache = context$abc_cache, filtered_meta = context$filtered_meta,
        gene_data = context$gene_data,
        selected_gene_data = context$selected_gene_data,
        download_base = context$download_base, data_type = context$data_type
    )
    context$selected_gene_data$.abc_gene_index <- NULL
    provenance <- .abc_output_provenance(
        context$filtered_meta, context$selected_gene_data,
        data_type = context$data_type, release = context$release
    )

    .abc_materialize(
        expression = expression, cell_metadata = context$filtered_meta,
        gene_data = context$selected_gene_data, data_type = context$data_type,
        return_type = return_type, coordinate_columns = coordinate_columns,
        provenance = provenance,
        spatial = isTRUE(attr(context$filtered_meta, "abc_spatial"))
    )
}

# Resolve the request from an approved plan, rejecting arguments the plan
# would otherwise silently override, or else from the individual arguments.
.abc_fetch_context <- function(
        plan, data, metadata, filters,
        gene_data, genes, data_type, offline,
        release, download_base) {
    if (is.null(plan)) {
        return(.abc_fetch_context_from_inputs(
            data = data, metadata = metadata, filters = filters,
            gene_data = gene_data, genes = genes, data_type = data_type,
            offline = offline, release = release, download_base = download_base
        ))
    }

    .abc_reject_plan_overrides(
        data = data, metadata = metadata, filters = filters,
        gene_data = gene_data, genes = genes, data_type = data_type,
        offline = offline, release = release, download_base = download_base
    )
    .abc_fetch_context_from_plan(plan)
}

#' Plan an ABC Atlas expression fetch without downloading expression files
#'
#' Applies the same cell filtering, gene selection and expression-file routing
#' as [fetch_data()], but does not download or read expression matrices.
#'
#' @param data Optional \code{abc_data} object returned by \code{load_data()}.
#'   This is the preferred interface. When supplied, its metadata, gene data,
#'   cache directory, and release are reused automatically.
#' @param download_base Directory where ABC Atlas data is downloaded and cached.
#' @param metadata Cell metadata, typically from [load_data()].
#' @param filters Named list of column = value filters.
#' @param gene_data Gene metadata, typically from [load_data()].
#' @param genes Character vector of genes to include. Each entry is matched
#'   first against \code{gene_identifier}, then against the published source
#'   identifier, then against \code{gene_symbol}; a symbol shared by several
#'   genes is reported as ambiguous. If \code{NULL}, all genes are used.
#' @param data_type Published expression representation (`"raw"` or `"log2"`).
#'   If `NULL`, the dataset route's declared default is used.
#' @param offline If `TRUE`, use only the locally cached manifest.
#' @param release ABC Atlas release. Defaults to the package-validated release.
#' @param confirm If `TRUE` (the default in interactive sessions) and a
#'   download is required, the plan is printed and the user is asked to
#'   approve it. In non-interactive sessions a plan that requires a download
#'   is returned unapproved.
#' @return An `abc_fetch_plan` object describing the requested output and
#'   expression files required to produce it.
#' @examples
#' example_dir <- system.file(
#'     "extdata", "ex",
#'     package = "ABCatlas"
#' )
#' example_data <- readRDS(
#'     file.path(example_dir, "example-data.rds")
#' )
#'
#' plan <- plan_fetch(
#'     download_base = example_dir,
#'     metadata = example_data$cells,
#'     gene_data = example_data$genes,
#'     genes = c("G1", "G4"),
#'     offline = TRUE,
#'     release = "20250101"
#' )
#' plan
#' @export
plan_fetch <- function(
        data = NULL, metadata = NULL, filters = list(),
        gene_data = NULL, genes = NULL, data_type = NULL,
        offline = FALSE, release = NULL,
        download_base = NULL,
        confirm = interactive()) {
    context <- .abc_fetch_context_from_inputs(
        data = data, metadata = metadata, filters = filters,
        gene_data = gene_data, genes = genes, data_type = data_type,
        offline = offline, release = release, download_base = download_base
    )
    files <- .abc_expression_files(
        context$abc_cache, context$filtered_meta,
        data_type = context$data_type
    )
    files$additional_download <- ifelse(
        files$cached, 0, pmax(files$size - files$partial_bytes, 0)
    )
    plan <- .abc_new_fetch_plan(context, files, offline)

    if (isTRUE(confirm) && plan$total_download > 0 && !isTRUE(offline) &&
        isTRUE(plan$storage_sufficient)) {
        print(plan)
        plan$approved <- .abc_confirm_download()
    }

    plan
}

# Assemble the abc_fetch_plan object: download totals, the disk-space check,
# and the resolved request that fetch_data(plan = ) reuses. A plan that needs
# no download is approved automatically.
.abc_new_fetch_plan <- function(context, files, offline) {
    total_download <- sum(files$additional_download)
    available_space <- .abc_available_space(context$download_base)

    structure(
        list(
            data_type = context$data_type,
            n_cells = nrow(context$filtered_meta),
            n_genes = nrow(context$selected_gene_data),
            files = files,
            total_size = sum(files$size),
            total_partial = sum(files$partial_bytes),
            total_download = total_download,
            available_space = available_space,
            storage_required = total_download,
            storage_sufficient = is.na(available_space) ||
                available_space >= total_download,
            release = context$abc_cache$version,
            approved = if (total_download == 0) TRUE else NA,
            .metadata = context$filtered_meta,
            .gene_data = context$gene_data,
            .selected_gene_data = context$selected_gene_data,
            .download_base = context$download_base,
            .offline = offline
        ),
        class = "abc_fetch_plan"
    )
}

#' @export
print.abc_fetch_plan <- function(x, ...) {
    cat(
        "<abc_fetch_plan>\n",
        "  release:        ", x$release, "\n",
        "  representation: ", x$data_type, "\n",
        "  requested:      ", x$n_genes, " genes x ", x$n_cells, " cells\n",
        "\nExpression files:\n",
        sep = ""
    )
    for (i in seq_len(nrow(x$files))) {
        .abc_print_plan_file(x$files[i, ])
    }

    if (is.na(x$available_space)) {
        space <- "unable to determine"
        check <- "unable to determine"
    } else {
        space <- .abc_format_size(x$available_space)
        check <- if (x$storage_sufficient) "PASS" else "INSUFFICIENT"
    }

    cat(
        "\nTotals:\n",
        "  selected-file size:   ", .abc_format_size(x$total_size), "\n",
        "  resumable partial:    ", .abc_format_size(x$total_partial), "\n",
        "  additional download:  ", .abc_format_size(x$total_download), "\n",
        "  available disk space: ", space, "\n",
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

.abc_validate_gene_resolution_inputs <- function(
        gene_data,
        genes
) {
    if (!"gene_identifier" %in% names(gene_data)) {
        stop(
            "'gene_identifier' column not found in gene_data.",
            call. = FALSE
        )
    }

    ids <- as.character(
        gene_data$gene_identifier
    )

    if (
        anyNA(ids) ||
            any(!nzchar(ids)) ||
            anyDuplicated(ids)
    ) {
        stop(
            "'gene_identifier' must be a non-missing unique key in gene_data.",
            call. = FALSE
        )
    }

    if (is.null(genes)) {
        return(
            list(
                ids = ids,
                query = NULL
            )
        )
    }

    query <- as.character(genes)

    if (
        !length(query) ||
            anyNA(query) ||
            any(!nzchar(query)) ||
            anyDuplicated(query)
    ) {
        stop(
            "genes must be a nonempty vector without duplicate requests.",
            call. = FALSE
        )
    }

    list(
        ids = ids,
        query = query
    )
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

.abc_resolve_genes <- function(
        gene_data,
        genes = NULL
) {
    inputs <- .abc_validate_gene_resolution_inputs(
        gene_data = gene_data,
        genes = genes
    )

    if (is.null(inputs$query)) {
        out <- gene_data
        out$.abc_gene_index <- seq_len(nrow(out))
        rownames(out) <- NULL
        return(out)
    }

    source_ids <- if (
        "source_gene_identifier" %in% names(gene_data)
    ) {
        as.character(
            gene_data$source_gene_identifier
        )
    } else {
        NULL
    }

    symbols <- if (
        "gene_symbol" %in% names(gene_data)
    ) {
        as.character(
            gene_data$gene_symbol
        )
    } else {
        NULL
    }

    selected <- vapply(
        inputs$query,
        .abc_resolve_gene_request,
        integer(1L),
        ids = inputs$ids,
        source_ids = source_ids,
        symbols = symbols,
        USE.NAMES = FALSE
    )

    .abc_finalize_gene_resolution(
        gene_data = gene_data,
        selected = selected
    )
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
            " absent after filtering: ",
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

.abc_select_cell_ids <- function(
        metadata,
        cell_ids
) {
    if (is.null(cell_ids)) {
        return(metadata)
    }

    labels <- .abc_validate_cell_id_selection(
        metadata = metadata,
        cell_ids = cell_ids
    )

    .abc_match_cell_ids(
        metadata = metadata,
        cell_ids = cell_ids,
        labels = labels
    )
}

.abc_apply_cell_filters <- function(
        metadata,
        filters
) {
    if (
        !is.list(filters) ||
            is.null(names(filters)) ||
            any(!nzchar(names(filters)))
    ) {
        if (length(filters)) {
            stop(
                "filters must be a named list of ",
                "metadata column = allowed values.",
                call. = FALSE
            )
        }
    }

    for (filter_col in names(filters)) {
        filter_val <- filters[[filter_col]]

        if (!filter_col %in% names(metadata)) {
            stop(
                "Column '",
                filter_col,
                "' not found in metadata.",
                call. = FALSE
            )
        }

        metadata <- metadata[
            metadata[[filter_col]] %in% filter_val, ,
            drop = FALSE
        ]
    }

    metadata
}

.abc_validate_cell_sample <- function(
        metadata,
        n_cells,
        seed
) {
    if (is.null(n_cells)) {
        return(NULL)
    }

    if (
        length(n_cells) != 1L ||
            is.na(n_cells) ||
            n_cells < 0 ||
            n_cells != as.integer(n_cells)
    ) {
        stop(
            "n_cells must be a single non-negative integer.",
            call. = FALSE
        )
    }

    n_cells <- as.integer(n_cells)

    if (n_cells > nrow(metadata)) {
        stop(
            "n_cells=",
            n_cells,
            " exceeds the ",
            nrow(metadata),
            " cells available after filtering.",
            call. = FALSE
        )
    }

    if (
        n_cells < nrow(metadata) &&
            (length(seed) != 1L || is.na(seed))
    ) {
        stop(
            "seed must be a single non-missing value.",
            call. = FALSE
        )
    }

    n_cells
}

.abc_sample_cell_rows <- function(
        metadata,
        n_cells,
        seed
) {
    keep <- withr::with_seed(
        seed,
        sample.int(
            nrow(metadata),
            size = n_cells,
            replace = FALSE
        )
    )

    metadata[
        keep, ,
        drop = FALSE
    ]
}

.abc_sample_cells <- function(
        metadata,
        n_cells,
        seed
) {
    n_cells <- .abc_validate_cell_sample(
        metadata = metadata,
        n_cells = n_cells,
        seed = seed
    )

    if (
        is.null(n_cells) ||
            n_cells == nrow(metadata)
    ) {
        return(metadata)
    }

    .abc_sample_cell_rows(
        metadata = metadata,
        n_cells = n_cells,
        seed = seed
    )
}

.abc_select_cells <- function(
        metadata, filters = list(), cell_ids = NULL,
        n_cells = NULL, seed = 1L) {
    out <- metadata

    out <- .abc_apply_cell_filters(
        metadata = out,
        filters = filters
    )

    out <- .abc_select_cell_ids(
        metadata = out,
        cell_ids = cell_ids
    )

    out <- .abc_sample_cells(
        metadata = out,
        n_cells = n_cells,
        seed = seed
    )

    if (!nrow(out)) {
        stop(
            "No cells left after selection. ",
            "Check filters, cell_ids or n_cells.",
            call. = FALSE
        )
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

# Apply cell filters and resolve the requested genes for a fetch.
.abc_prepare_fetch <- function(metadata, filters, gene_data, genes = NULL) {
    filtered_meta <- .abc_select_cells(
        metadata,
        filters = filters
    )

    resolved_genes <- .abc_resolve_genes(
        gene_data,
        genes
    )

    list(
        metadata = filtered_meta,
        gene_data = resolved_genes
    )
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
        cached = file.exists(attrs$local_path),
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
