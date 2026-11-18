test_that("abc_matrix prints a short summary", {
    example_dir <- system.file("extdata", "ex", package = "ABCatlas")
    example_data <- readRDS(file.path(example_dir, "example-data.rds"))
    x <- suppressMessages(fetch_data(
        download_base = example_dir, metadata = example_data$cells,
        gene_data = example_data$genes, genes = c("G1", "G4"),
        offline = TRUE, release = "20250101"
    ))
    out <- capture.output(res <- print(x))
    expect_identical(out[1], "<abc_matrix>")
    expect_true(any(grepl("2 genes x 6 cells", out, fixed = TRUE)))
    expect_lt(length(out), 10)
    expect_identical(res, x)
})
