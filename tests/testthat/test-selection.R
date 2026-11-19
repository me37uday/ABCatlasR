test_that("fetch_data returns the cells in the order of cell_ids", {
    fx <- make_fixture_cache()
    x <- suppressMessages(fetch_data(
        fixture_data(fx),
        cell_ids = c("c4", "c1"),
        genes = c("G1", "G4"),
        offline = TRUE
    ))
    expect_identical(colnames(x$expression), c("c4", "c1"))
    expect_identical(x$cell_metadata$cell_label, c("c4", "c1"))
    expect_equal(
        unname(as.matrix(x$expression)),
        unname(fx$truth[c("G1", "G4"), c("c4", "c1")])
    )
})

test_that("fetch_data with genes = NULL returns every gene", {
    fx <- make_fixture_cache()
    x <- suppressMessages(fetch_data(
        fixture_data(fx),
        cell_ids = fx$cells$cell_label,
        offline = TRUE
    ))
    expect_identical(rownames(x$expression), fx$genes$gene_identifier)
    expect_equal(unname(as.matrix(x$expression)), unname(fx$truth))
})

test_that("fetch_data does not need plan_fetch", {
    fx <- make_fixture_cache()
    local_mocked_bindings(
        plan_fetch = function(...) stop("plan_fetch was called", call. = FALSE)
    )
    x <- suppressMessages(fetch_data(
        fixture_data(fx),
        cell_ids = "c2",
        offline = TRUE
    ))
    expect_identical(colnames(x$expression), "c2")
})

test_that("the three functions take the documented arguments", {
    expect_identical(
        names(formals(load_data)),
        c("dataset", "download_base", "offline", "release")
    )
    expect_identical(
        names(formals(plan_fetch)),
        c("data", "cell_ids", "data_type", "offline")
    )
    expect_identical(
        names(formals(fetch_data)),
        c("data", "cell_ids", "genes", "data_type", "return_type", "offline")
    )
})
