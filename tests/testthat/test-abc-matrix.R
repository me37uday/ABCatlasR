test_that("abc_matrix prints a short summary", {
    abc <- example_abc_data()
    x <- suppressMessages(fetch_data(
        abc,
        cell_ids = abc$cell_metadata$cell_label,
        genes = c("G1", "G4"),
        offline = TRUE
    ))
    out <- capture.output(res <- print(x))
    expect_identical(out[1], "<abc_matrix>")
    expect_true(any(grepl("2 genes x 6 cells", out, fixed = TRUE)))
    expect_lt(length(out), 10)
    expect_identical(res, x)
})
