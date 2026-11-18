#' @export
print.abc_matrix <- function(x, ...) {
    p <- x$provenance
    cat(
        "<abc_matrix>\n",
        "  dataset:        ", p$dataset, "\n",
        "  release:        ", p$release, "\n",
        "  representation: ", x$data_type, "\n",
        "  expression:     ", nrow(x$expression), " genes x ",
        ncol(x$expression), " cells (sparse)\n",
        "  elements:       ", paste(names(x), collapse = ", "), "\n",
        sep = ""
    )
    invisible(x)
}
