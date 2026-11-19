# Rebuild the offline example release shipped in inst/extdata/ex.
#
# The data are synthetic: six cells, four genes and two expression files
# (TEST-A and TEST-B, each as raw and log2), laid out like an ABC Atlas
# download directory with a release manifest, metadata CSV files and H5AD
# expression matrices. They are not derived from ABC Atlas data. The man
# page examples of abc_cache_connect(), plan_fetch() and fetch_data(), and
# some unit tests, read them with offline = TRUE.
#
# make_fixture_cache() in tests/testthat/helper-fixtures.R builds the same
# files, except that it writes the taxonomy to "TEST-taxonomy" instead of
# "tax".
#
# Usage, from the package root (the output directory must be new or empty):
#   Rscript inst/scripts/make-example-data.R <output directory>
#
# Requires the Matrix, rhdf5 and jsonlite packages.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
    stop("Usage: Rscript make-example-data.R <output directory>")
}
out <- args[1]
if (dir.exists(out) && length(list.files(out))) {
    stop(out, " is not empty.")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
version <- "20250101"

# Write a minimal h5ad holding `counts` (cells x genes) as CSR, in the
# on-disk layout used by the ABC Atlas files (anndata >= 0.8).
write_h5ad <- function(path, counts, cell_labels, gene_ids, gene_symbols) {
    csr <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "RsparseMatrix")
    rhdf5::h5createFile(path)
    fid <- rhdf5::H5Fopen(path)
    on.exit(rhdf5::H5Fclose(fid))

    set_attrs <- function(group, attrs) {
        gid <- rhdf5::H5Gopen(fid, group)
        on.exit(rhdf5::H5Gclose(gid))
        for (n in names(attrs)) rhdf5::h5writeAttribute(attrs[[n]], gid, n)
    }

    rhdf5::h5createGroup(fid, "X")
    rhdf5::h5write(as.numeric(csr@x), fid, "X/data")
    rhdf5::h5write(as.integer(csr@j), fid, "X/indices")
    rhdf5::h5write(as.integer(csr@p), fid, "X/indptr")
    set_attrs("X", list(
        `encoding-type` = "csr_matrix", `encoding-version` = "0.1.0",
        shape = as.integer(dim(counts))
    ))

    rhdf5::h5createGroup(fid, "obs")
    rhdf5::h5write(cell_labels, fid, "obs/cell_label")
    set_attrs("obs", list(
        `_index` = "cell_label", `encoding-type` = "dataframe"
    ))

    rhdf5::h5createGroup(fid, "var")
    rhdf5::h5write(gene_ids, fid, "var/gene_identifier")
    rhdf5::h5write(gene_symbols, fid, "var/gene_symbol")
    set_attrs("var", list(
        `_index` = "gene_identifier", `encoding-type` = "dataframe"
    ))
    invisible(path)
}

# Record a written file in the release manifest listing.
listing <- list()
add_file <- function(directory, kind, name, rel, sub = NULL) {
    path <- file.path(out, rel)
    entry <- list(
        version = version, relative_path = rel,
        url = paste0("https://example.invalid/", rel),
        size = file.size(path), file_hash = unname(tools::md5sum(path))
    )
    files <- list(files = stats::setNames(list(entry), tools::file_ext(rel)))
    if (is.null(sub)) {
        listing[[directory]][[kind]][[name]] <<- files
    } else {
        listing[[directory]][[kind]][[name]][[sub]] <<- files
    }
}

write_csv <- function(df, directory, name) {
    rel <- file.path("metadata", directory, version, paste0(name, ".csv"))
    dir.create(
        file.path(out, dirname(rel)),
        recursive = TRUE, showWarnings = FALSE
    )
    utils::write.csv(df, file.path(out, rel), row.names = FALSE, na = "")
    add_file(directory, "metadata", name, rel)
}

# Cell and gene metadata. The symbol "Aaa" is shared by two genes, and the
# trailing space in "alpha " is kept, so that both cases can be exercised.
cells <- data.frame(
    cell_label = c("c1", "c2", "c3", "c4", "c5", "c6"),
    feature_matrix_label = rep(c("TEST-A", "TEST-B"), each = 3),
    dataset_label = "TEST-10X",
    donor_sex = c("M", "F", "M", "F", NA, "M"),
    region_of_interest_label = c("R1", "R2", "R1", "R2", "R1", "R1"),
    cluster_alias = c(1L, 2L, 1L, 2L, 2L, 1L)
)
write_csv(cells, "TEST-10X", "cell_metadata")

genes <- data.frame(
    gene_identifier = c("G1", "G2", "G3", "G4"),
    gene_symbol = c("Aaa", "Bbb", "Aaa", "Ddd"),
    name = c("alpha ", "beta", "gamma", "delta")
)
write_csv(genes, "TEST-10X", "gene")

write_csv(
    data.frame(
        region_of_interest_label = c("R1", "R2"),
        color_hex_triplet = c("#111111", "#222222")
    ),
    "TEST-10X", "region_of_interest_structure_map"
)

# Taxonomy: two clusters, each with a cluster and a supercluster term.
write_csv(
    data.frame(
        cluster_alias = c(1L, 1L, 2L, 2L),
        cluster_annotation_term_set_name = c(
            "cluster", "supercluster", "cluster", "supercluster"
        ),
        cluster_annotation_term_name = c("cl1", "sup1", "cl2", NA),
        color_hex_triplet = c("#aa0000", "#bb0000", "#cc0000", "#dd0000")
    ),
    "tax", "cluster_to_cluster_annotation_membership"
)
write_csv(
    data.frame(label = c("L1", "L2"), name = c("supercluster", "cluster")),
    "tax", "cluster_annotation_term_set"
)

# Expression (genes x cells), split across two files by cell. Within each
# file the cells are stored in reverse order, so that the on-disk order
# differs from the metadata order.
truth <- matrix(
    c(
        1, 0, 3, 0,
        0, 2, 0, 0,
        5, 0, 0, 1,
        0, 0, 7, 0,
        0, 4, 0, 0,
        0, 0, 0, 0
    ),
    nrow = 4,
    dimnames = list(genes$gene_identifier, cells$cell_label)
)
truth_log2 <- log2(truth + 1)

for (f in c("TEST-A", "TEST-B")) {
    in_file <- rev(cells$cell_label[cells$feature_matrix_label == f])
    for (data_type in c("raw", "log2")) {
        values <- if (data_type == "raw") truth else truth_log2
        rel <- file.path(
            "expression_matrices", "TEST-10X", version,
            paste0(f, "-", data_type, ".h5ad")
        )
        dir.create(
            file.path(out, dirname(rel)),
            recursive = TRUE, showWarnings = FALSE
        )
        write_h5ad(
            file.path(out, rel), t(values[, in_file]), in_file,
            genes$gene_identifier, genes$gene_symbol
        )
        add_file("TEST-10X", "expression_matrices", f, rel, sub = data_type)
    }
}

# Release manifest.
manifest <- list(
    version = version, resource_uri = "s3://example.invalid/",
    directory_listing = list(), file_listing = listing
)
dir.create(file.path(out, "releases", version), recursive = TRUE)
jsonlite::write_json(
    manifest, file.path(out, "releases", version, "manifest.json"),
    auto_unbox = TRUE
)

# example-data.rds: an abc_data object for these files, in the form
# load_data() returns, used by the man page examples. Its download_base is
# NULL; the examples point it at the installed example directory. The
# missing donor_sex is stored as an empty string, as in the CSV file.
cells$donor_sex[is.na(cells$donor_sex)] <- ""
rownames(cells) <- cells$cell_label
attr(cells, "abc_dataset") <- "TEST"
attr(cells, "abc_default_data_type") <- "raw"
attr(cells, "abc_available_data_types") <- c("raw", "log2")
attr(cells, "abc_expression_directories") <- "TEST-10X"
attr(cells, "abc_modality") <- "single_cell_transcriptomics"
attr(cells, "abc_spatial") <- FALSE
example <- list(
    cell_metadata = cells,
    gene_data = genes,
    unique_values = lapply(
        cells[c("donor_sex", "region_of_interest_label")], unique
    ),
    dataset = "TEST",
    download_base = NULL,
    release = version
)
class(example) <- c("abc_data", "list")
saveRDS(example, file.path(out, "example-data.rds"))
