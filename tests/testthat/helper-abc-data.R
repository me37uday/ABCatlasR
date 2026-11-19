# The fixture as the abc_data object load_data() would return for it.
fixture_data <- function(fx, cells = fx$cells, genes = fx$genes) {
    rownames(cells) <- cells$cell_label
    structure(
        list(
            cell_metadata = cells,
            gene_data = genes,
            unique_values = list(),
            dataset = "TEST",
            download_base = fx$dir,
            release = fx$version
        ),
        class = c("abc_data", "list")
    )
}

# The example data shipped in inst/extdata/ex, pointed at its files.
example_abc_data <- function() {
    example_dir <- system.file("extdata", "ex", package = "ABCatlasR")
    abc <- readRDS(file.path(example_dir, "example-data.rds"))
    abc$download_base <- example_dir
    abc
}
