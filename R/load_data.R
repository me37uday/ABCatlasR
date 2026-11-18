# =============================================================================
# Unified metadata loader.
#
# Datasets differ only in declarative details (which directory holds the cell
# metadata, which auxiliary tables to join and on what keys, whether there is a
# taxonomy/region step, and the columns of interest). Those details live in the
# .abc_datasets registry below; .build_dataset() turns a spec into the assembled
# result, and load_data() is the single user-facing entry point.
# =============================================================================

# ---- Dataset registry -------------------------------------------------------
# Each entry is a declarative spec. Fields:
#   cell             list(dir, file = "cell_metadata",
#                         col_classes = <named chr>)
#   taxonomy         taxonomy directory (triggers the cluster reshape + join)
#   roi              list(dir, file) for the region-of-interest colour join
#   joins            ordered list of list(dir, file, key, suffix) left-joins
#   gene             list(dir, file = "gene")
#   cols_of_interest columns to summarise into `unique_values`
.abc_datasets <- list(
    WHB = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "WHB-10Xv3",
        modality = "single_nucleus_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "WHB-10Xv3", col_classes = c(cell_label = "character")
        ),
        taxonomy = "WHB-taxonomy",
        roi = list(
            dir = "WHB-10Xv3", file = "region_of_interest_structure_map"
        ),
        gene = list(dir = "WHB-10Xv3"),
        cols_of_interest = c(
            "feature_matrix_label", "brain_section_label",
            "region_of_interest_label", "anatomical_division_label",
            "donor_sex",
            "subcluster", "cluster", "supercluster", "neurotransmitter"
        )
    ),
    WMB = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = c(
            "WMB-10Xv2", "WMB-10Xv3",
            "WMB-10XMulti"
        ),
        modality = "single_cell_transcriptomics",
        spatial = FALSE,
        cell = list(dir = "WMB-10X", col_classes = c(cell_label = "character")),
        taxonomy = "WMB-taxonomy",
        gene = list(dir = "WMB-10X"),
        cols_of_interest = c(
            "feature_matrix_label", "library_method",
            "region_of_interest_acronym", "donor_sex", "class", "subclass",
            "supertype", "cluster", "neurotransmitter"
        )
    ),
    ConsensusMouseAIBS = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = c(
            "WMB-10Xv2",
            "WMB-10Xv3",
            "WMB-10XMulti"
        ),
        modality = "single_cell_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "Consensus-WMB-AIBS-10X",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "Consensus-WMB-integrated-taxonomy",
        gene = list(dir = "WMB-10X"),
        cols_of_interest = character()
    ),
    ConsensusMouseMacosko = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "Consensus-WMB-Macosko-10X",
        modality = "single_nucleus_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "Consensus-WMB-Macosko-10X",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "Consensus-WMB-integrated-taxonomy",
        gene = list(dir = "Consensus-WMB-Macosko-10X"),
        cols_of_interest = character()
    ),
    DevelopingMouse = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "Developing-Mouse-Vis-Cortex-10X",
        modality = "single_cell_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "Developing-Mouse-Vis-Cortex-10X",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "Developing-Mouse-Vis-Cortex-taxonomy",
        gene = list(dir = "Developing-Mouse-Vis-Cortex-10X"),
        cols_of_interest = c("synchronized_age")
    ),
    AgingMouse = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "Zeng-Aging-Mouse-10Xv3",
        modality = "single_cell_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "Zeng-Aging-Mouse-10Xv3",
            col_classes = c(
                cell_label = "character",
                wmb_cluster_alias = "integer"
            )
        ),
        # Match the validated Python-reticulate annotation path: native Aging
        # annotations define cluster identity and age bias, while the separate
        # mapping directory contributes the optional WMB reference mapping.
        joins = list(
            list(
                dir = "Zeng-Aging-Mouse-10Xv3",
                file = "cell_cluster_annotations",
                key = "cell_label"
            ),
            list(
                dir = "Zeng-Aging-Mouse-WMB-taxonomy",
                file = "cell_cross_mapping_annotations",
                key = "cell_label"
            )
        ),
        # No gene table of its own: the expression matrix uses the WMB 10X gene
        # set (same 32285 genes, same order).
        gene = list(dir = "WMB-10X"),
        cols_of_interest = c(
            "anatomical_division_label",
            "donor_sex",
            "donor_age",
            "cluster_name",
            "cluster_age_bias",
            "neurotransmitter_combined_label",
            "wmb_class_name",
            "wmb_subclass_name",
            "wmb_supertype_name",
            "wmb_cluster_name"
        )
    ),
    HMBA = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "HMBA-10xMultiome-BG-Aligned",
        modality = "single_nucleus_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "HMBA-10xMultiome-BG-Aligned",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "HMBA-BG-taxonomy-CCN20250428",
        joins = list(
            list(
                dir = "HMBA-10xMultiome-BG-Aligned",
                file = "donor",
                key = "donor_label"
            ),
            list(
                dir = "HMBA-10xMultiome-BG-Aligned",
                file = "library",
                key = "library_label",
                column_suffix = "_library_table"
            )
        ),
        gene = list(dir = "HMBA-10xMultiome-BG-Aligned"),
        cols_of_interest = c(
            "species_scientific_name",
            "species_common_name",
            "donor_sex",
            "donor_age",
            "region_of_interest_name",
            "anatomical_division_label",
            "Neighborhood",
            "Class",
            "Subclass",
            "Group",
            "Cluster"
        )
    ),
    HMBA_Human = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "HMBA-10xMultiome-BG",
        modality = "single_nucleus_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "HMBA-10xMultiome-BG",
            col_classes = c(cell_label = "character")
        ),
        subset = list(
            dir = "HMBA-10xMultiome-BG",
            file = "donor",
            key = "donor_label",
            column = "species_genus",
            value = "Human"
        ),
        taxonomy = "HMBA-BG-taxonomy-CCN20250428",
        gene = list(
            dir = "HMBA-10xMultiome-BG",
            file = "human_gene"
        ),
        cols_of_interest = c(
            "species_genus",
            "species_scientific_name",
            "species_common_name",
            "donor_sex",
            "donor_age",
            "region_of_interest_name",
            "anatomical_division_label"
        )
    ),
    HMBA_Macaque = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "HMBA-10xMultiome-BG",
        modality = "single_nucleus_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "HMBA-10xMultiome-BG",
            col_classes = c(cell_label = "character")
        ),
        subset = list(
            dir = "HMBA-10xMultiome-BG",
            file = "donor",
            key = "donor_label",
            column = "species_genus",
            value = "Macaque"
        ),
        taxonomy = "HMBA-BG-taxonomy-CCN20250428",
        gene = list(
            dir = "HMBA-10xMultiome-BG",
            file = "macaque_gene"
        ),
        cols_of_interest = c(
            "species_genus",
            "species_scientific_name",
            "species_common_name",
            "donor_sex",
            "donor_age",
            "region_of_interest_name",
            "anatomical_division_label"
        )
    ),
    HMBA_Marmoset = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "HMBA-10xMultiome-BG",
        modality = "single_nucleus_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "HMBA-10xMultiome-BG",
            col_classes = c(cell_label = "character")
        ),
        subset = list(
            dir = "HMBA-10xMultiome-BG",
            file = "donor",
            key = "donor_label",
            column = "species_genus",
            value = "Marmoset"
        ),
        taxonomy = "HMBA-BG-taxonomy-CCN20250428",
        gene = list(
            dir = "HMBA-10xMultiome-BG",
            file = "marmoset_gene"
        ),
        cols_of_interest = c(
            "species_genus",
            "species_scientific_name",
            "species_common_name",
            "donor_sex",
            "donor_age",
            "region_of_interest_name",
            "anatomical_division_label"
        )
    ),
    HMBA_PatchSeq = list(
        default_data_type = "log2",
        available_data_types = "log2",
        expression_directories = "HMBA-Macaque-PatchSeq",
        modality = "patch_seq_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "HMBA-Macaque-PatchSeq",
            col_classes = c(cell_label = "character")
        ),
        auto_joins = TRUE,
        taxonomy = "HMBA-BG-taxonomy-CCN20250428",
        gene = list(dir = "HMBA-Macaque-PatchSeq"),
        cols_of_interest = c(
            "library_technique",
            "barcoded_cell_sample_technique",
            "region_of_interest_label",
            "donor_age",
            "donor_sex",
            "species_genus",
            "species_scientific_name",
            "Neighborhood",
            "Class",
            "Subclass",
            "Group",
            "Cluster"
        )
    ),
    HMBA_spatial_Human = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "HMBA-MERSCOPE-H22.30.001-BG",
        modality = "spatial_transcriptomics",
        spatial = TRUE,
        cell = list(
            dir = "HMBA-MERSCOPE-H22.30.001-BG",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "HMBA-BG-taxonomy-CCN20250428",
        gene = list(dir = "HMBA-MERSCOPE-H22.30.001-BG"),
        cols_of_interest = c(
            "brain_section_label",
            "x_experiment",
            "y_experiment",
            "x_slab_mm",
            "y_slab_mm",
            "parcellation_term_symbol",
            "Neighborhood",
            "Class",
            "Subclass",
            "Group",
            "Cluster"
        )
    ),
    HMBA_spatial_Macaque = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "HMBA-MERSCOPE-QM23.50.001-BG",
        modality = "spatial_transcriptomics",
        spatial = TRUE,
        cell = list(
            dir = "HMBA-MERSCOPE-QM23.50.001-BG",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "HMBA-BG-taxonomy-CCN20250428",
        gene = list(dir = "HMBA-MERSCOPE-QM23.50.001-BG"),
        cols_of_interest = c(
            "brain_section_label",
            "x_experiment",
            "y_experiment",
            "x_slab_mm",
            "y_slab_mm",
            "parcellation_term_symbol",
            "Neighborhood",
            "Class",
            "Subclass",
            "Group",
            "Cluster"
        )
    ),
    HMBA_spatial_Marmoset = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "HMBA-Xenium-CJ23.56.004-BG",
        modality = "spatial_transcriptomics",
        spatial = TRUE,
        cell = list(
            dir = "HMBA-Xenium-CJ23.56.004-BG",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "HMBA-BG-taxonomy-CCN20250428",
        gene = list(dir = "HMBA-Xenium-CJ23.56.004-BG"),
        cols_of_interest = c(
            "brain_section_label",
            "x_experiment",
            "y_experiment",
            "x_rotated_um",
            "y_rotated_um",
            "x_slab_mm",
            "y_slab_mm",
            "parcellation_term_symbol",
            "Neighborhood",
            "Class",
            "Subclass",
            "Group",
            "Cluster"
        )
    ),
    Zhuang1 = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "Zhuang-ABCA-1",
        modality = "spatial_transcriptomics",
        spatial = TRUE,
        cell = list(
            dir = "Zhuang-ABCA-1",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "WMB-taxonomy",
        gene = list(dir = "Zhuang-ABCA-1"),
        cols_of_interest = c(
            "brain_section_label",
            "donor_sex",
            "donor_genotype",
            "x", "y", "z",
            "subclass",
            "supertype",
            "cluster",
            "neurotransmitter",
            "subclass_confidence_score",
            "cluster_confidence_score",
            "high_quality_transfer"
        )
    ),
    Zhuang2 = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "Zhuang-ABCA-2",
        modality = "spatial_transcriptomics",
        spatial = TRUE,
        cell = list(
            dir = "Zhuang-ABCA-2",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "WMB-taxonomy",
        gene = list(dir = "Zhuang-ABCA-2"),
        cols_of_interest = c(
            "brain_section_label",
            "donor_sex",
            "donor_genotype",
            "x", "y", "z",
            "subclass",
            "supertype",
            "cluster",
            "neurotransmitter",
            "subclass_confidence_score",
            "cluster_confidence_score",
            "high_quality_transfer"
        )
    ),
    Zhuang3 = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "Zhuang-ABCA-3",
        modality = "spatial_transcriptomics",
        spatial = TRUE,
        cell = list(
            dir = "Zhuang-ABCA-3",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "WMB-taxonomy",
        gene = list(dir = "Zhuang-ABCA-3"),
        cols_of_interest = c(
            "brain_section_label",
            "donor_sex",
            "donor_genotype",
            "x", "y", "z",
            "subclass",
            "supertype",
            "cluster",
            "neurotransmitter",
            "subclass_confidence_score",
            "cluster_confidence_score",
            "high_quality_transfer"
        )
    ),
    Zhuang4 = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "Zhuang-ABCA-4",
        modality = "spatial_transcriptomics",
        spatial = TRUE,
        cell = list(
            dir = "Zhuang-ABCA-4",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "WMB-taxonomy",
        gene = list(dir = "Zhuang-ABCA-4"),
        cols_of_interest = c(
            "brain_section_label",
            "donor_sex",
            "donor_genotype",
            "x", "y", "z",
            "subclass",
            "supertype",
            "cluster",
            "neurotransmitter",
            "subclass_confidence_score",
            "cluster_confidence_score",
            "high_quality_transfer"
        )
    ),
    MERFISH = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "MERFISH-C57BL6J-638850",
        modality = "spatial_transcriptomics",
        spatial = TRUE,
        cell = list(
            dir = "MERFISH-C57BL6J-638850",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "WMB-taxonomy",
        gene = list(dir = "MERFISH-C57BL6J-638850"),
        cols_of_interest = c(
            "brain_section_label",
            "donor_sex",
            "donor_genotype",
            "x",
            "y",
            "z",
            "class",
            "subclass",
            "supertype",
            "cluster",
            "neurotransmitter"
        )
    ),
    MERFISH_imputed = list(
        default_data_type = "log2",
        available_data_types = "log2",
        expression_directories = "MERFISH-C57BL6J-638850-imputed",
        expression_matrix_map = c(
            "C57BL6J-638850" = "C57BL6J-638850-imputed"
        ),
        modality = "imputed_spatial_transcriptomics",
        spatial = TRUE,
        cell = list(
            dir = "MERFISH-C57BL6J-638850",
            col_classes = c(cell_label = "character")
        ),
        taxonomy = "WMB-taxonomy",
        gene = list(dir = "WMB-10X"),
        cols_of_interest = c(
            "brain_section_label",
            "donor_sex",
            "donor_genotype",
            "x",
            "y",
            "z",
            "class",
            "subclass",
            "supertype",
            "cluster",
            "neurotransmitter"
        )
    ),
    SEAAD_Multiregion = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "SEA-AD-Multiregion-10X",
        modality = "single_nucleus_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "SEA-AD-Multiregion-10X",
            col_classes = c(cell_label = "character")
        ),
        auto_joins = TRUE,
        taxonomy = "SEA-AD-Multiregion-taxonomy",
        gene = list(dir = "SEA-AD-Multiregion-10X"),
        cols_of_interest = c(
            "feature_matrix_label",
            "region_of_interest_label",
            "region_of_interest_name",
            "Brain Region",
            "species_scientific_name",
            "species_genus",
            "donor_sex",
            "donor_age",
            "donor_race",
            "Overall AD neuropathological Change",
            "Braak",
            "CERAD score",
            "Cognitive Status",
            "APOE Genotype"
        )
    ),
    SEAAD_CaH = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "SEA-AD-CaH-10X",
        modality = "single_nucleus_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "SEA-AD-CaH-10X",
            col_classes = c(cell_label = "character")
        ),
        auto_joins = TRUE,
        taxonomy = "SEA-AD-CaH-taxonomy",
        gene = list(dir = "SEA-AD-CaH-10X"),
        cols_of_interest = c(
            "feature_matrix_label",
            "region_of_interest_label",
            "region_of_interest_name",
            "Brain Region",
            "species_scientific_name",
            "species_genus",
            "donor_sex",
            "donor_age",
            "donor_race",
            "Overall AD neuropathological Change",
            "Braak",
            "CERAD score",
            "Cognitive Status",
            "APOE Genotype"
        )
    ),
    PMDBS = list(
        default_data_type = "raw",
        available_data_types = c("raw", "log2"),
        expression_directories = "ASAP-PMDBS-10X",
        modality = "single_nucleus_transcriptomics",
        spatial = FALSE,
        cell = list(
            dir = "ASAP-PMDBS-10X",
            col_classes = c(cell_label = "character")
        ),
        joins = list(
            list(dir = "ASAP-PMDBS-10X", file = "sample", key = "sample_label"),
            list(dir = "ASAP-PMDBS-10X", file = "donor", key = "donor_label")
        ),
        # PMDBS expression matrices use the native 36601-gene human gene
        # universe.
        gene = list(dir = "ASAP-PMDBS-10X"),
        cols_of_interest = c(
            "source_dataset_label", "region_of_interest_label",
            "donor_race", "donor_sex", "primary_diagnosis", "age_at_death",
            "apoe4_status", "braak_stage", "cerad_score", "cognitive_status"
        )
    )
)

#' Load assembled cell metadata for an ABC Atlas dataset
#'
#' Single entry point for every supported dataset. The dataset is selected by
#' name; what differs between datasets (directories, joins, taxonomy, columns of
#' interest) is described declaratively in an internal registry rather than in
#' per-dataset functions.
#'
#' @param dataset Dataset route name. Run [list_datasets()] to see every
#'   supported route (for example \code{"WHB"}, \code{"AgingMouse"} or
#'   \code{"MERFISH"}). Unambiguous partial names are accepted.
#' @param download_base Directory where ABC Atlas data is downloaded and cached.
#'   If \code{NULL} (default), it is resolved by \code{\link{abc_download_dir}}.
#' @param offline If \code{TRUE}, use only locally cached data and never contact
#'   S3 (see \code{\link{abc_cache_connect}}).
#' @param release Optional ABC Atlas release version in \code{YYYYMMDD} form.
#'   When supplied, metadata are loaded from exactly that release.
#' @param filters Named list of metadata column = allowed values. Filters are
#'   applied before exact cell selection and sampling.
#' @param cell_ids Optional character vector of exact cell labels. Requested
#'   order is preserved and missing cells are an error.
#' @param n_cells Optional number of cells to sample after filters and
#'   \code{cell_ids} selection.
#' @param seed Random seed used when \code{n_cells} samples the selected cells.
#' @return An \code{abc_data} object containing \code{cell_metadata},
#'   \code{gene_data}, \code{unique_values}, and dataset/cache provenance.
#' @examples
#' # Cache the downloaded metadata in a temporary folder.
#' abc <- load_data(
#'     "HMBA_PatchSeq",
#'     download_base = file.path(tempdir(), "ABCatlas")
#' )
#' abc
#' dim(abc$cell_metadata)
#' @export
load_data <- function(
        dataset, download_base = NULL, offline = FALSE,
        release = .abc_default_release, filters = list(),
        cell_ids = NULL, n_cells = NULL, seed = 1L) {
    datasets <- names(.abc_datasets)
    if (missing(dataset)) {
        stop("Specify a dataset, one of: ", paste(datasets, collapse = ", "))
    }
    dataset <- match.arg(dataset, datasets)
    spec <- .abc_datasets[[dataset]]

    download_base <- abc_download_dir(download_base)
    abc_cache <- abc_cache_connect(
        download_base,
        offline = offline, release = release
    )

    out <- .build_dataset(spec, abc_cache, dataset)
    out$cell_metadata <- .abc_select_cells(
        out$cell_metadata,
        filters = filters, cell_ids = cell_ids,
        n_cells = n_cells, seed = seed
    )

    # Recompute filter choices from the selected population while retaining
    # the route's curated columns of interest.
    cols <- intersect(
        spec$cols_of_interest %||% character(), names(out$cell_metadata)
    )
    out$unique_values <- lapply(out$cell_metadata[cols], unique)

    out$dataset <- dataset
    out$download_base <- download_base
    out$release <- abc_cache$version
    class(out) <- c("abc_data", "list")
    out
}

#' @export
print.abc_data <- function(x, ...) {
    metadata <- x$cell_metadata
    gene_data <- x$gene_data

    default_type <- attr(
        metadata,
        "abc_default_data_type"
    ) %||% "raw"

    available <- attr(
        metadata,
        "abc_available_data_types"
    ) %||% c("raw", "log2")

    modality <- attr(
        metadata,
        "abc_modality"
    ) %||% "transcriptomics"

    spatial <- isTRUE(
        attr(metadata, "abc_spatial")
    )

    cat(
        "<abc_data>\n",
        "  dataset:        ", x$dataset, "\n",
        "  release:        ", x$release, "\n",
        "  cells:          ", nrow(metadata), "\n",
        "  genes:          ", nrow(gene_data), "\n",
        "  modality:       ", modality, "\n",
        "  spatial:        ", if (spatial) "yes" else "no", "\n",
        "  expression:     ", paste(available, collapse = ", "),
        " (default: ", default_type, ")\n",
        "  expression data: not loaded\n",
        sep = ""
    )

    invisible(x)
}

# ---- Generic builder --------------------------------------------------------

# Left many-to-one metadata join that preserves cell count and row order.
# Annotation tables must have a non-missing unique join key; otherwise joining
# them could silently duplicate cells.
.abc_safe_join <- function(
        left, right, key, suffix = "annotation",
        column_suffix = paste0("__", suffix)) {
    if (!key %in% names(left) || !key %in% names(right)) {
        stop(suffix, ": join key ", key, " missing.", call. = FALSE)
    }

    left[[key]] <- as.character(left[[key]])
    right[[key]] <- as.character(right[[key]])
    right <- unique(right)

    if (anyNA(right[[key]]) || anyDuplicated(right[[key]])) {
        stop(
            suffix, ": right-hand ", key,
            " is not a non-missing unique key; join would multiply cells.",
            call. = FALSE
        )
    }

    n_before <- nrow(left)
    left$.abc_join_order <- seq_len(n_before)

    overlap <- setdiff(intersect(names(left), names(right)), key)
    if (length(overlap)) {
        replacement <- paste0(overlap, column_suffix)

        if (any(replacement %in% union(names(left), names(right)))) {
            stop(suffix, ": ambiguous annotation column suffix.", call. = FALSE)
        }

        names(right)[match(overlap, names(right))] <- replacement
    }

    out <- dplyr::left_join(left, right, by = key)
    out <- out[order(out$.abc_join_order), , drop = FALSE]
    out$.abc_join_order <- NULL

    if (nrow(out) != n_before) {
        stop("Join changed the number of cells.", call. = FALSE)
    }

    out
}

# Make published gene identifiers unique without dropping or aggregating
# features. Duplicate source IDs retain their original identity and occurrence
# so expression readers can reconcile them one-to-one against the H5AD axis.
.abc_normalize_gene_metadata <- function(gene, context = "gene metadata") {
    if (!"gene_identifier" %in% names(gene)) {
        stop(context, ": required column 'gene_identifier' missing.",
            call. = FALSE
        )
    }

    out <- gene
    original <- as.character(out$gene_identifier)
    if (anyNA(original) || any(!nzchar(original))) {
        stop(context, ": missing gene_identifier.", call. = FALSE)
    }

    counts <- table(original, useNA = "no")
    if (!any(counts > 1L)) {
        out$gene_identifier <- original
        return(out)
    }

    if ("source_gene_identifier" %in% names(out)) {
        stop(context, ": source_gene_identifier already exists; ",
            "cannot safely disambiguate duplicates.",
            call. = FALSE
        )
    }

    ids <- .abc_disambiguate_gene_ids(original, counts)
    out$source_gene_identifier <- original
    out$source_gene_occurrence <- ids$occurrence
    out$source_gene_group_size <- as.integer(counts[original])
    out$gene_identifier <- ids$unique_ids

    if (anyNA(out$gene_identifier) || any(!nzchar(out$gene_identifier)) ||
        anyDuplicated(out$gene_identifier)) {
        stop(context,
            ": duplicate gene identifiers remain after normalization; ",
            "no silent aggregation allowed.",
            call. = FALSE
        )
    }

    out
}

# Give each member of a duplicated identifier group a unique "<id>.dup<k>"
# name, skipping suffixes that already exist as identifiers. `occurrence` is
# the true ordinal within the original group; it deliberately differs from k
# when an existing identifier forced a suffix to be skipped.
.abc_disambiguate_gene_ids <- function(original, counts) {
    used <- unique(original)
    seen <- integer()
    unique_ids <- character(length(original))
    occurrence <- integer(length(original))

    for (i in seq_along(original)) {
        value <- original[[i]]
        if (counts[[value]] <= 1L) {
            unique_ids[[i]] <- value
            occurrence[[i]] <- 1L
            next
        }

        index <- (if (value %in% names(seen)) seen[[value]] else 0L) + 1L
        candidate <- paste0(value, ".dup", index)
        while (candidate %in% used) {
            index <- index + 1L
            candidate <- paste0(value, ".dup", index)
        }

        seen[[value]] <- index
        used <- c(used, candidate)
        unique_ids[[i]] <- candidate
        occurrence[[i]] <- sum(original[seq_len(i)] == value)
    }

    list(unique_ids = unique_ids, occurrence = occurrence)
}

# Turn a dataset spec into the assembled list(cell_metadata, gene_data,
# unique_values). All dataset-specific behaviour comes from `spec`.
.abc_dataset_subset_keys <- function(subset_metadata, subset_spec, dataset) {
    missing <- setdiff(
        c(subset_spec$key, subset_spec$column), names(subset_metadata)
    )
    if (length(missing)) {
        stop(
            dataset, ": subset metadata is missing required column(s): ",
            paste(missing, collapse = ", "), ".",
            call. = FALSE
        )
    }

    values <- subset_metadata[[subset_spec$column]]
    matched <- !is.na(values) & as.character(values) == subset_spec$value
    keep_keys <- unique(as.character(subset_metadata[matched, subset_spec$key]))
    keep_keys <- keep_keys[!is.na(keep_keys) & nzchar(keep_keys)]

    if (!length(keep_keys)) {
        stop(
            dataset, ": no ", subset_spec$value, " cells found from ",
            subset_spec$file, ".",
            call. = FALSE
        )
    }

    keep_keys
}

.abc_apply_dataset_subset <- function(
        cell_metadata,
        subset_spec,
        abc_cache,
        dataset
) {
    subset_metadata <- abc_cache$get_metadata_dataframe(
        directory = subset_spec$dir,
        file_name = subset_spec$file
    )

    keep_keys <- .abc_dataset_subset_keys(
        subset_metadata = subset_metadata,
        subset_spec = subset_spec,
        dataset = dataset
    )

    cell_keys <- as.character(
        cell_metadata[[subset_spec$key]]
    )

    cell_metadata <- cell_metadata[
        !is.na(cell_keys) &
            cell_keys %in% keep_keys, ,
        drop = FALSE
    ]

    .abc_safe_join(
        cell_metadata,
        subset_metadata,
        subset_spec$key,
        subset_spec$file
    )
}

.abc_apply_standard_metadata_joins <- function(
        cell_metadata,
        cell_dir,
        abc_cache
) {
    standard_joins <- list(
        list(file = "library", key = "library_label"),
        list(file = "sample", key = "sample_label"),
        list(file = "donor", key = "donor_label"),
        list(file = "disease", key = "donor_label"),
        list(
            file = "specimen_metadata",
            key = "brain_section_label"
        )
    )

    for (join in standard_joins) {
        if (
            join$key %in% names(cell_metadata) &&
                .abc_metadata_exists(
                    abc_cache,
                    cell_dir,
                    join$file
                )
        ) {
            auxiliary <- abc_cache$get_metadata_dataframe(
                directory = cell_dir,
                file_name = join$file
            )

            cell_metadata <- .abc_safe_join(
                cell_metadata,
                auxiliary,
                join$key,
                join$file
            )
        }
    }

    cell_metadata
}

.abc_apply_taxonomy_membership <- function(
        cell_metadata, cell_dir,
        taxonomy_dir, abc_cache) {
    membership_dir <- NULL
    for (candidate in c(cell_dir, taxonomy_dir)) {
        if (.abc_metadata_exists(
            abc_cache, candidate,
            "cell_to_cluster_membership"
        )) {
            membership_dir <- candidate
            break
        }
    }

    if (!is.null(membership_dir)) {
        membership <- abc_cache$get_metadata_dataframe(
            directory = membership_dir, file_name = "cell_to_cluster_membership"
        )
        if ("cluster_alias" %in% names(cell_metadata)) {
            names(cell_metadata)[names(cell_metadata) == "cluster_alias"] <-
                "cluster_alias_original"
            cell_metadata$cluster_alias_original <-
                as.character(cell_metadata$cluster_alias_original)
        }
        cell_metadata <- .abc_safe_join(
            cell_metadata, membership, "cell_label", "membership"
        )
    }

    if (!"cluster_alias" %in% names(cell_metadata)) {
        stop(cell_dir, ": no cluster_alias or applicable membership table.",
            call. = FALSE
        )
    }

    cell_metadata$cluster_alias <- as.character(cell_metadata$cluster_alias)
    cell_metadata
}

.abc_join_taxonomy_annotations <- function(
        cell_metadata,
        taxonomy_dir,
        abc_cache
) {
    taxonomy <- .build_taxonomy(
        abc_cache,
        taxonomy_dir
    )

    cell_metadata <- .abc_safe_join(
        cell_metadata,
        taxonomy$cluster_details,
        "cluster_alias",
        "taxonomy"
    )

    if (!is.null(taxonomy$cluster_colors)) {
        cell_metadata <- .abc_safe_join(
            cell_metadata,
            taxonomy$cluster_colors,
            "cluster_alias",
            "taxonomy_color",
            column_suffix = "_color"
        )
    }

    cell_metadata
}

.abc_apply_dataset_taxonomy <- function(
        cell_metadata,
        cell_dir,
        taxonomy_dir,
        abc_cache
) {
    cell_metadata <- .abc_apply_taxonomy_membership(
        cell_metadata = cell_metadata,
        cell_dir = cell_dir,
        taxonomy_dir = taxonomy_dir,
        abc_cache = abc_cache
    )

    .abc_join_taxonomy_annotations(
        cell_metadata = cell_metadata,
        taxonomy_dir = taxonomy_dir,
        abc_cache = abc_cache
    )
}

.abc_apply_spatial_metadata <- function(
        cell_metadata,
        cell_dir,
        abc_cache
) {
    spatial_files <- c(
        "slab_plane_coordinates",
        "cell_anatomical_annotations"
    )

    for (file in spatial_files) {
        if (.abc_metadata_exists(
            abc_cache,
            cell_dir,
            file
        )) {
            spatial_metadata <- abc_cache$get_metadata_dataframe(
                directory = cell_dir,
                file_name = file
            )

            cell_metadata <- .abc_safe_join(
                cell_metadata,
                spatial_metadata,
                "cell_label",
                file
            )
        }
    }

    cell_metadata
}

.abc_apply_roi_colors <- function(
        cell_metadata,
        roi_spec,
        abc_cache
) {
    roi <- abc_cache$get_metadata_dataframe(
        directory = roi_spec$dir,
        file_name = roi_spec$file
    )

    roi$region_of_interest_label <- make.unique(
        as.character(roi$region_of_interest_label)
    )

    names(roi)[
        names(roi) == "color_hex_triplet"
    ] <- "region_of_interest_color"

    roi <- roi[
        ,
        c(
            "region_of_interest_label",
            "region_of_interest_color"
        )
    ]

    dplyr::left_join(
        cell_metadata,
        roi,
        by = "region_of_interest_label"
    )
}

.abc_apply_configured_metadata_joins <- function(
        cell_metadata,
        joins,
        abc_cache
) {
    for (join in joins %||% list()) {
        auxiliary <- abc_cache$get_metadata_dataframe(
            directory = join$dir,
            file_name = join$file
        )

        cell_metadata <- .abc_safe_join(
            cell_metadata,
            auxiliary,
            join$key,
            join$file,
            column_suffix = join$column_suffix %||%
                paste0("__", join$file)
        )
    }

    cell_metadata
}

.abc_load_dataset_genes <- function(
        gene_spec,
        abc_cache,
        dataset
) {
    gene <- abc_cache$get_metadata_dataframe(
        directory = gene_spec$dir,
        file_name = gene_spec$file %||% "gene"
    )

    gene <- .abc_normalize_gene_metadata(
        gene,
        context = gene_spec$dir
    )

    rownames(gene) <- gene$gene_identifier

    message(
        sprintf(
            "%s: %s genes loaded.",
            dataset,
            format(nrow(gene), big.mark = ",")
        )
    )

    gene
}

.abc_finalize_dataset <- function(
        cell_metadata,
        gene_data,
        spec,
        dataset
) {
    cols <- intersect(
        spec$cols_of_interest,
        colnames(cell_metadata)
    )

    unique_values <- lapply(
        cell_metadata[cols],
        unique
    )

    attr(cell_metadata, "abc_dataset") <- dataset
    attr(
        cell_metadata,
        "abc_default_data_type"
    ) <- spec$default_data_type
    attr(
        cell_metadata,
        "abc_available_data_types"
    ) <- spec$available_data_types
    attr(
        cell_metadata,
        "abc_expression_directories"
    ) <- spec$expression_directories
    attr(
        cell_metadata,
        "abc_expression_matrix_map"
    ) <- spec$expression_matrix_map
    attr(
        cell_metadata,
        "abc_modality"
    ) <- spec$modality
    attr(
        cell_metadata,
        "abc_spatial"
    ) <- spec$spatial

    list(
        cell_metadata = cell_metadata,
        gene_data = gene_data,
        unique_values = unique_values
    )
}

.abc_enrich_dataset_cells <- function(
        cell_metadata,
        spec,
        abc_cache
) {
    if (!is.null(spec$taxonomy)) {
        cell_metadata <- .abc_apply_dataset_taxonomy(
            cell_metadata = cell_metadata,
            cell_dir = spec$cell$dir,
            taxonomy_dir = spec$taxonomy,
            abc_cache = abc_cache
        )
    }

    if (isTRUE(spec$spatial)) {
        cell_metadata <- .abc_apply_spatial_metadata(
            cell_metadata = cell_metadata,
            cell_dir = spec$cell$dir,
            abc_cache = abc_cache
        )
    }

    if (!is.null(spec$roi)) {
        cell_metadata <- .abc_apply_roi_colors(
            cell_metadata = cell_metadata,
            roi_spec = spec$roi,
            abc_cache = abc_cache
        )
    }

    cell_metadata <- .abc_apply_configured_metadata_joins(
        cell_metadata = cell_metadata,
        joins = spec$joins,
        abc_cache = abc_cache
    )

    message(
        "Cell annotations loaded."
    )

    cell_metadata
}

.abc_assemble_dataset_cells <- function(
        spec,
        abc_cache,
        dataset
) {
    cell_metadata <- abc_cache$get_metadata_dataframe(
        directory = spec$cell$dir,
        file_name = spec$cell$file %||% "cell_metadata",
        col_classes = spec$cell$col_classes
    )

    if (!is.null(spec$subset)) {
        cell_metadata <- .abc_apply_dataset_subset(
            cell_metadata = cell_metadata,
            subset_spec = spec$subset,
            abc_cache = abc_cache,
            dataset = dataset
        )
    }

    if (isTRUE(spec$auto_joins)) {
        cell_metadata <- .abc_apply_standard_metadata_joins(
            cell_metadata = cell_metadata,
            cell_dir = spec$cell$dir,
            abc_cache = abc_cache
        )
    }

    message(
        sprintf(
            "%s: %s cells loaded.",
            dataset,
            format(nrow(cell_metadata), big.mark = ",")
        )
    )

    .abc_enrich_dataset_cells(
        cell_metadata = cell_metadata,
        spec = spec,
        abc_cache = abc_cache
    )
}

.build_dataset <- function(spec, abc_cache, dataset = "") {
    cell_metadata <- .abc_assemble_dataset_cells(
        spec = spec,
        abc_cache = abc_cache,
        dataset = dataset
    )

    gene_data <- .abc_load_dataset_genes(
        gene_spec = spec$gene,
        abc_cache = abc_cache,
        dataset = dataset
    )

    .abc_finalize_dataset(
        cell_metadata = cell_metadata,
        gene_data = gene_data,
        spec = spec,
        dataset = dataset
    )
}

# Build the wide cluster_details / cluster_colors lookup tables from a taxonomy
# directory's membership + term-set tables, keyed by cluster_alias. Shared by
# every dataset that declares a `taxonomy` (currently WHB and WMB).
.abc_validate_taxonomy_alias <- function(
        x,
        source,
        taxonomy_dir
) {
    if (!"cluster_alias" %in% names(x)) {
        stop(
            taxonomy_dir,
            ": ",
            source,
            " is missing cluster_alias.",
            call. = FALSE
        )
    }

    x$cluster_alias <- as.character(
        x$cluster_alias
    )

    if (
        anyNA(x$cluster_alias) ||
            anyDuplicated(x$cluster_alias)
    ) {
        stop(
            taxonomy_dir,
            ": ",
            source,
            " does not have a non-missing unique cluster_alias.",
            call. = FALSE
        )
    }

    x
}

.abc_add_taxonomy_term_names <- function(
        membership,
        abc_cache,
        taxonomy_dir
) {
    if ("cluster_annotation_term_name" %in% names(membership)) {
        return(membership)
    }

    if (!.abc_metadata_exists(
        abc_cache,
        taxonomy_dir,
        "cluster_annotation_term"
    )) {
        stop(
            taxonomy_dir,
            ": term names absent and no term table available.",
            call. = FALSE
        )
    }

    term <- abc_cache$get_metadata_dataframe(
        directory = taxonomy_dir,
        file_name = "cluster_annotation_term"
    )

    names(term)[
        names(term) == "label"
    ] <- "cluster_annotation_term_label"

    names(term)[
        names(term) == "name"
    ] <- "cluster_annotation_term_name"

    .abc_safe_join(
        membership,
        term[
            ,
            c(
                "cluster_annotation_term_label",
                "cluster_annotation_term_name"
            ),
            drop = FALSE
        ],
        "cluster_annotation_term_label",
        "term"
    )
}

.abc_add_taxonomy_term_set_names <- function(
        membership,
        abc_cache,
        taxonomy_dir
) {
    if ("cluster_annotation_term_set_name" %in% names(membership)) {
        return(membership)
    }

    if (!.abc_metadata_exists(
        abc_cache,
        taxonomy_dir,
        "cluster_annotation_term_set"
    )) {
        stop(
            taxonomy_dir,
            ": term-set names absent and no term-set table available.",
            call. = FALSE
        )
    }

    term_sets <- abc_cache$get_metadata_dataframe(
        directory = taxonomy_dir,
        file_name = "cluster_annotation_term_set"
    )

    names(term_sets)[
        names(term_sets) == "label"
    ] <- "cluster_annotation_term_set_label"

    names(term_sets)[
        names(term_sets) == "name"
    ] <- "cluster_annotation_term_set_name"

    .abc_safe_join(
        membership,
        term_sets[
            ,
            c(
                "cluster_annotation_term_set_label",
                "cluster_annotation_term_set_name"
            ),
            drop = FALSE
        ],
        "cluster_annotation_term_set_label",
        "termset"
    )
}

.abc_enrich_taxonomy_membership <- function(
        membership,
        abc_cache,
        taxonomy_dir
) {
    membership <- .abc_add_taxonomy_term_names(
        membership = membership,
        abc_cache = abc_cache,
        taxonomy_dir = taxonomy_dir
    )

    .abc_add_taxonomy_term_set_names(
        membership = membership,
        abc_cache = abc_cache,
        taxonomy_dir = taxonomy_dir
    )
}

.abc_validate_taxonomy_membership_keys <- function(
        membership,
        taxonomy_dir
) {
    required <- c(
        "cluster_alias",
        "cluster_annotation_term_set_name",
        "cluster_annotation_term_name"
    )

    missing <- setdiff(
        required,
        names(membership)
    )

    if (length(missing)) {
        stop(
            taxonomy_dir,
            ": missing taxonomy field(s): ",
            paste(missing, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    membership$cluster_alias <- as.character(
        membership$cluster_alias
    )

    if (
        anyNA(membership$cluster_alias) ||
            anyNA(membership$cluster_annotation_term_set_name)
    ) {
        stop(
            taxonomy_dir,
            ": missing taxonomy keys.",
            call. = FALSE
        )
    }

    membership
}

.abc_validate_taxonomy_membership_conflicts <- function(
        membership,
        taxonomy_dir
) {
    key <- paste(
        membership$cluster_alias,
        membership$cluster_annotation_term_set_name,
        sep = "\r"
    )

    conflicts <- tapply(
        membership$cluster_annotation_term_name,
        key,
        function(x) length(unique(x)) > 1L
    )

    if (any(conflicts)) {
        stop(
            taxonomy_dir,
            ": conflicting cluster annotations.",
            call. = FALSE
        )
    }

    membership
}

.abc_validate_taxonomy_membership <- function(
        membership,
        taxonomy_dir
) {
    membership <- .abc_validate_taxonomy_membership_keys(
        membership = membership,
        taxonomy_dir = taxonomy_dir
    )

    .abc_validate_taxonomy_membership_conflicts(
        membership = membership,
        taxonomy_dir = taxonomy_dir
    )
}

.abc_load_taxonomy_membership <- function(
        abc_cache,
        taxonomy_dir
) {
    if (!.abc_metadata_exists(
        abc_cache,
        taxonomy_dir,
        "cluster_to_cluster_annotation_membership"
    )) {
        return(NULL)
    }

    membership <- abc_cache$get_metadata_dataframe(
        directory = taxonomy_dir,
        file_name = "cluster_to_cluster_annotation_membership"
    )

    membership <- .abc_enrich_taxonomy_membership(
        membership = membership,
        abc_cache = abc_cache,
        taxonomy_dir = taxonomy_dir
    )

    .abc_validate_taxonomy_membership(
        membership = membership,
        taxonomy_dir = taxonomy_dir
    )
}

.abc_load_pivoted_taxonomy_details <- function(
        abc_cache,
        taxonomy_dir
) {
    .abc_validate_taxonomy_alias(
        abc_cache$get_metadata_dataframe(
            directory = taxonomy_dir,
            file_name = "cluster_to_cluster_annotation_membership_pivoted"
        ),
        "pivoted taxonomy",
        taxonomy_dir
    )
}

.abc_reshape_taxonomy_details <- function(
        membership,
        taxonomy_dir
) {
    if (is.null(membership)) {
        stop(
            taxonomy_dir,
            ": taxonomy membership file missing.",
            call. = FALSE
        )
    }

    detail <- unique(
        membership[
            ,
            c(
                "cluster_alias",
                "cluster_annotation_term_set_name",
                "cluster_annotation_term_name"
            )
        ]
    )

    cluster_details <- stats::reshape(
        detail,
        idvar = "cluster_alias",
        timevar = "cluster_annotation_term_set_name",
        direction = "wide"
    )

    names(cluster_details) <- sub(
        "^cluster_annotation_term_name\\.",
        "",
        names(cluster_details)
    )

    .abc_validate_taxonomy_alias(
        cluster_details,
        "taxonomy annotations",
        taxonomy_dir
    )
}

.abc_get_taxonomy_details <- function(
        abc_cache,
        taxonomy_dir,
        membership
) {
    if (.abc_metadata_exists(
        abc_cache,
        taxonomy_dir,
        "cluster_to_cluster_annotation_membership_pivoted"
    )) {
        return(
            .abc_load_pivoted_taxonomy_details(
                abc_cache = abc_cache,
                taxonomy_dir = taxonomy_dir
            )
        )
    }

    .abc_reshape_taxonomy_details(
        membership = membership,
        taxonomy_dir = taxonomy_dir
    )
}

.abc_fill_taxonomy_details <- function(cluster_details) {
    detail_cols <- setdiff(
        names(cluster_details),
        "cluster_alias"
    )

    for (column in detail_cols) {
        cluster_details[[column]][
            is.na(cluster_details[[column]])
        ] <- "Other"
    }

    cluster_details
}

.abc_build_taxonomy_details <- function(
        abc_cache,
        taxonomy_dir,
        membership
) {
    cluster_details <- .abc_get_taxonomy_details(
        abc_cache = abc_cache,
        taxonomy_dir = taxonomy_dir,
        membership = membership
    )

    .abc_fill_taxonomy_details(
        cluster_details
    )
}

.abc_validate_taxonomy_color_conflicts <- function(
        color,
        taxonomy_dir
) {
    color_key <- paste(
        color$cluster_alias,
        color$cluster_annotation_term_set_name,
        sep = "\r"
    )

    color_conflicts <- tapply(
        color$color_hex_triplet,
        color_key,
        function(x) length(unique(x[!is.na(x)])) > 1L
    )

    if (any(color_conflicts)) {
        stop(
            taxonomy_dir,
            ": conflicting taxonomy colours.",
            call. = FALSE
        )
    }

    color
}

.abc_reshape_taxonomy_colors <- function(
        color,
        taxonomy_dir
) {
    cluster_colors <- stats::reshape(
        color,
        idvar = "cluster_alias",
        timevar = "cluster_annotation_term_set_name",
        direction = "wide"
    )

    names(cluster_colors) <- sub(
        "^color_hex_triplet\\.",
        "",
        names(cluster_colors)
    )

    .abc_validate_taxonomy_alias(
        cluster_colors,
        "taxonomy colours",
        taxonomy_dir
    )
}

.abc_build_taxonomy_colors <- function(
        membership,
        taxonomy_dir
) {
    if (
        is.null(membership) ||
            !"color_hex_triplet" %in% names(membership)
    ) {
        return(NULL)
    }

    color <- unique(
        membership[
            ,
            c(
                "cluster_alias",
                "cluster_annotation_term_set_name",
                "color_hex_triplet"
            )
        ]
    )

    color <- .abc_validate_taxonomy_color_conflicts(
        color = color,
        taxonomy_dir = taxonomy_dir
    )

    .abc_reshape_taxonomy_colors(
        color = color,
        taxonomy_dir = taxonomy_dir
    )
}

.build_taxonomy <- function(abc_cache, taxonomy_dir) {
    membership <- .abc_load_taxonomy_membership(
        abc_cache = abc_cache,
        taxonomy_dir = taxonomy_dir
    )

    cluster_details <- .abc_build_taxonomy_details(
        abc_cache = abc_cache,
        taxonomy_dir = taxonomy_dir,
        membership = membership
    )

    cluster_colors <- .abc_build_taxonomy_colors(
        membership = membership,
        taxonomy_dir = taxonomy_dir
    )

    list(
        cluster_details = cluster_details,
        cluster_colors = cluster_colors
    )
}
