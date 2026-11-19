plan_args <- function(fx) {
    list(
        download_base = fx$dir, metadata = fx$cells, gene_data = fx$genes,
        genes = c("G1", "G4"), offline = TRUE, release = fx$version
    )
}

test_that("plan_fetch samples cells reproducibly and selects exact cells", {
    fx <- make_fixture_cache()
    args <- plan_args(fx)
    p1 <- suppressMessages(do.call(plan_fetch, c(args, n_cells = 2, seed = 1)))
    p2 <- suppressMessages(do.call(plan_fetch, c(args, n_cells = 2, seed = 1)))
    expect_identical(p1$n_cells, 2L)
    expect_identical(rownames(p1$.metadata), rownames(p2$.metadata))

    p3 <- suppressMessages(
        do.call(plan_fetch, c(args, list(cell_ids = c("c4", "c1"))))
    )
    expect_identical(rownames(p3$.metadata), c("c4", "c1"))
})

test_that("fetch_data returns exactly the planned cells", {
    fx <- make_fixture_cache()
    args <- c(plan_args(fx), n_cells = 3, seed = 2, confirm = FALSE)
    plan <- suppressMessages(do.call(plan_fetch, args))
    x <- suppressMessages(fetch_data(plan))
    expect_identical(colnames(x$expression), rownames(plan$.metadata))
})

test_that("confirm = FALSE approves a plan that needs a download", {
    fx <- make_fixture_cache()
    cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
    unlink(.abc_file_attributes(cache, "TEST-10X", "TEST-A/raw")$local_path)
    plan <- suppressMessages(
        do.call(plan_fetch, c(plan_args(fx), confirm = FALSE))
    )
    expect_gt(plan$total_download, 0)
    expect_true(plan$approved)
})

test_that("fetch_data takes a plan and load_data no longer selects cells", {
    expect_error(fetch_data(list()), "abc_fetch_plan")
    expect_false(any(
        c("filters", "cell_ids", "n_cells", "seed") %in% names(formals(load_data))
    ))
})
