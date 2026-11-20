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

test_that("cell selection preserves ABC routing provenance", {
    metadata <- data.frame(
        cell_label = paste0("c", 1:5),
        group = c("A", "A", "B", "B", "B"),
        stringsAsFactors = FALSE
    )

    attr(metadata, "abc_dataset") <- "TEST"
    attr(metadata, "abc_default_data_type") <- "log2"
    attr(metadata, "abc_available_data_types") <- "log2"
    attr(metadata, "abc_expression_directories") <- "TEST-DIR"
    attr(metadata, "abc_modality") <- "test_transcriptomics"
    attr(metadata, "abc_spatial") <- FALSE

    out <- .abc_select_cells(
        metadata,
        cell_ids = c("c5", "c3")
    )

    expect_identical(
        attr(out, "abc_dataset"),
        "TEST"
    )
    expect_identical(
        attr(out, "abc_default_data_type"),
        "log2"
    )
    expect_identical(
        attr(out, "abc_available_data_types"),
        "log2"
    )
    expect_identical(
        attr(out, "abc_expression_directories"),
        "TEST-DIR"
    )
    expect_identical(
        attr(out, "abc_modality"),
        "test_transcriptomics"
    )
    expect_identical(
        attr(out, "abc_spatial"),
        FALSE
    )
})
