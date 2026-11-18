dense <- function(m) as.matrix(m)

test_that("expression file planning matches cell routing without downloading", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  files <- .abc_expression_files(
    cache,
    fx$cells,
    data_type = "raw"
  )

  expect_equal(
    files$feature_matrix_label,
    c("TEST-A", "TEST-B")
  )

  expect_equal(files$n_cells, c(3, 3))
  expect_true(all(files$data_type == "raw"))
  expect_true(all(files$size > 0))
  expect_true(all(files$cached))
  expect_equal(files$partial_bytes, c(0, 0))
})

test_that("expression file planning reports resumable partial bytes", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  attrs <- .abc_file_attributes(
    cache,
    "TEST-10X",
    "TEST-A/raw"
  )

  part <- paste0(attrs$local_path, ".part")
  writeBin(as.raw(1:10), part)

  files <- .abc_expression_files(
    cache,
    fx$cells,
    data_type = "raw"
  )

  expect_equal(
    files$partial_bytes[
      files$feature_matrix_label == "TEST-A"
    ],
    10
  )
})

test_that("expression file planning validates routing columns", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  bad <- fx$cells
  bad$feature_matrix_label <- NULL

  expect_error(
    .abc_expression_files(cache, bad),
    "feature_matrix_label"
  )
})

test_that("expression matrix matches the truth for all cells and genes", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  cells <- fx$cells
  rownames(cells) <- cells$cell_label

  m <- suppressMessages(
    .abc_expression_matrix(
      cache,
      cells,
      fx$genes,
      seq_len(nrow(fx$genes))
    )
  )
  expect_s4_class(m, "dgCMatrix")
  expect_equal(rownames(m), c("G1", "G2", "G3", "G4"))
  expect_equal(colnames(m), cells$cell_label)
  expect_equal(unname(dense(m)), unname(fx$truth))
})

test_that("gene and cell selection, order and block size are respected", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  cells <- fx$cells[c(5, 1, 6, 3), ]          # mixed files, not on-disk order
  rownames(cells) <- cells$cell_label

  m <- suppressMessages(
    .abc_expression_matrix(
      cache,
      cells,
      fx$genes,
      c(1L, 3L, 4L),
      block_size = 1L
    )
  )
  expect_equal(rownames(m), c("G1", "G3", "G4"))
  expect_equal(colnames(m), c("c5", "c1", "c6", "c3"))
  expect_equal(unname(dense(m)), unname(fx$truth[c("G1", "G3", "G4"), c("c5", "c1", "c6", "c3")]))
})

test_that("cells absent from their expression file give empty columns", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  cells <- rbind(fx$cells[1:2, ], transform(fx$cells[1, ], cell_label = "ghost"))
  rownames(cells) <- cells$cell_label

  m <- suppressMessages(
    .abc_expression_matrix(
      cache,
      cells,
      fx$genes,
      seq_len(nrow(fx$genes))
    )
  )
  expect_equal(colnames(m), c("c1", "c2", "ghost"))
  expect_equal(sum(m[, "ghost"]), 0)
  expect_equal(unname(dense(m[, 1:2])), unname(fx$truth[, 1:2]))
})

test_that("a gene table from another dataset is reported", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  cells <- fx$cells
  rownames(cells) <- cells$cell_label
  other <- transform(fx$genes, gene_identifier = paste0("X", gene_identifier))
  expect_error(
    suppressMessages(
      .abc_expression_matrix(
        cache,
        cells,
        other,
        c(1L, 3L)
      )
    ),
    "missing from"
  )
})

test_that("malformed identifiers in the file fall back to positional alignment", {
  # As in Zeng-Aging-Mouse-10Xv3: a printed pandas Series instead of the ID.
  garbled <- "gene_symbol\nAaa    G1\nAaa    G3\nName: gene_identifier, dtype: object"
  fx <- make_fixture_cache(var_ids = c("G1", "G2", garbled, "G4"))
  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  cells <- fx$cells
  rownames(cells) <- cells$cell_label

  warns <- capture_warnings(
    m <- suppressMessages(
      .abc_expression_matrix(
        cache,
        cells,
        fx$genes,
        c(1L, 3L)
      )
    )
  )
  expect_length(warns, 2)   # one per expression file
  expect_match(warns, "1 gene identifiers in TEST-[AB]-raw.h5ad are malformed")
  expect_equal(unname(dense(m)), unname(fx$truth[c("G1", "G3"), ]))
  # Genes whose identifiers are intact need no fallback.
  expect_no_warning(
    suppressMessages(
      .abc_expression_matrix(
        cache,
        cells,
        fx$genes,
        4L
      )
    )
  )
})

test_that("fetch_data filters cells and returns the default matrix bundle", {

  fx <- make_fixture_cache()

  obj <- suppressMessages(
    fetch_data(
      download_base = fx$dir,
      metadata = fx$cells,
      filters = list(donor_sex = "M"),
      gene_data = fx$genes,
      genes = c("G1", "G3", "G4"),
      offline = TRUE,
      release = fx$version
    )
  )

  expect_s3_class(
    obj,
    "abc_matrix"
  )

  expect_equal(
    colnames(obj$expression),
    c("c1", "c3", "c6")
  )

  expect_equal(
    obj$cell_metadata$donor_sex,
    c("M", "M", "M")
  )

  expect_equal(
    unname(dense(obj$expression)),
    unname(
      fx$truth[
        c("G1", "G3", "G4"),
        c("c1", "c3", "c6")
      ]
    )
  )
})

test_that("fetch_data validates filters", {
  fx <- make_fixture_cache()
  expect_error(fetch_data(download_base = fx$dir, metadata = fx$cells,
                          filters = list(nope = 1), gene_data = fx$genes,
                          offline = TRUE, release = fx$version),
               "not found in metadata")
  expect_error(fetch_data(download_base = fx$dir, metadata = fx$cells,
                          filters = list(donor_sex = "X"), gene_data = fx$genes,
                          offline = TRUE, release = fx$version),
               "No cells left")
})

test_that("fetch_data stops before download when disk space is insufficient", {
  fx <- make_fixture_cache()

  attrs <- .abc_file_attributes(
    abc_cache_connect(
      fx$dir,
      offline = TRUE,
      release = fx$version
    ),
    "TEST-10X",
    "TEST-A/raw"
  )

  unlink(attrs$local_path)

  local_mocked_bindings(
    .abc_available_space = function(path) 1,
    .abc_expression_matrix = function(...) {
      stop("expression matrix should not be read", call. = FALSE)
    }
  )

  expect_error(
    fetch_data(
      download_base = fx$dir,
      metadata = fx$cells,
      filters = list(donor_sex = "M"),
      gene_data = fx$genes,
      genes = c("G1", "G3", "G4"),
      offline = TRUE,
      release = fx$version
    ),
    "Insufficient disk space"
  )
})

test_that("plan_fetch reports the request without downloading expression files", {
  fx <- make_fixture_cache()

  plan <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    filters = list(donor_sex = "M"),
    gene_data = fx$genes,
    genes = c("G1", "G3", "G4"),
    offline = TRUE,
    release = fx$version
  )

  expect_s3_class(plan, "abc_fetch_plan")
  expect_equal(plan$data_type, "raw")
  expect_equal(plan$n_cells, 3)
  expect_equal(plan$n_genes, 3)
  expect_equal(plan$release, fx$version)
  expect_equal(plan$files$feature_matrix_label, c("TEST-A", "TEST-B"))
  expect_equal(plan$files$n_cells, c(2, 1))
  expect_equal(plan$total_download, 0)
  expect_equal(plan$storage_required, 0)
  expect_type(plan$available_space, "double")
  expect_length(plan$available_space, 1L)
  expect_type(plan$storage_sufficient, "logical")
  expect_length(plan$storage_sufficient, 1L)
})

test_that("approved fetch plan executes without a second confirmation", {
  fx <- make_fixture_cache()

  plan <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    gene_data = fx$genes,
    genes = c("G1", "G3", "G4"),
    offline = TRUE,
    release = fx$version,
    confirm = FALSE
  )

  expect_true(plan$approved)

  local_mocked_bindings(
    .abc_confirm_download = function(...) {
      stop("fetch_data should not ask for confirmation", call. = FALSE)
    }
  )

  out <- fetch_data(
    plan = plan,
    return_type = "matrix"
  )

  expect_s3_class(out, "abc_matrix")
  expect_equal(dim(out$expression), c(3L, 6L))
  expect_identical(
    rownames(out$expression),
    c("G1", "G3", "G4")
  )
  expect_identical(
    colnames(out$expression),
    paste0("c", 1:6)
  )
})

test_that("unapproved fetch plan is rejected", {
  fx <- make_fixture_cache()

  plan <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    gene_data = fx$genes,
    genes = "G1",
    offline = TRUE,
    release = fx$version,
    confirm = FALSE
  )

  plan$approved <- FALSE

  expect_error(
    fetch_data(plan = plan),
    "fetch plan was not approved"
  )
})

test_that("plan_fetch reports remaining bytes without downloading", {
  fx <- make_fixture_cache()

  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  attrs <- .abc_file_attributes(cache, "TEST-10X", "TEST-A/raw")

  unlink(attrs$local_path)
  part <- paste0(attrs$local_path, ".part")
  writeBin(as.raw(1:10), part)

  size_before <- file.info(part)$size

  plan <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    gene_data = fx$genes,
    genes = "G1",
    offline = TRUE,
    release = fx$version
  )

  row <- plan$files[plan$files$feature_matrix_label == "TEST-A", ]

  expect_false(row$cached)
  expect_equal(row$partial_bytes, 10)
  expect_equal(row$additional_download, row$size - 10)
  expect_equal(file.info(part)$size, size_before)
  expect_false(file.exists(attrs$local_path))
})

test_that("abc_fetch_plan prints a readable preflight summary", {
  fx <- make_fixture_cache()

  plan <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    filters = list(donor_sex = "M"),
    gene_data = fx$genes,
    genes = c("G1", "G3", "G4"),
    offline = TRUE,
    release = fx$version
  )

  expect_output(print(plan), "<abc_fetch_plan>")
  expect_output(print(plan), "3 genes x 3 cells")
  expect_output(print(plan), "TEST-10X/TEST-A/raw")
  expect_output(print(plan), "additional download")
  expect_output(print(plan), "available disk space")
  expect_output(print(plan), "storage check")
  expect_output(print(plan), "No expression files were downloaded")
})

test_that("data type resolution uses the dataset default", {
  fx <- make_fixture_cache()

  expect_equal(
    .abc_resolve_data_type(NULL, fx$cells),
    "raw"
  )
})

test_that("data type resolution accepts published representations", {
  fx <- make_fixture_cache()

  expect_equal(.abc_resolve_data_type("raw", fx$cells), "raw")
  expect_equal(.abc_resolve_data_type("log2", fx$cells), "log2")
})

test_that("data type resolution rejects invalid representations", {
  fx <- make_fixture_cache()

  expect_error(
    .abc_resolve_data_type("banana", fx$cells),
    "data_type must be one of"
  )
})

test_that("data type resolution rejects unavailable representations", {
  fx <- make_fixture_cache()
  attr(fx$cells, "abc_default_data_type") <- "log2"
  attr(fx$cells, "abc_available_data_types") <- "log2"

  expect_error(
    .abc_resolve_data_type("raw", fx$cells),
    "does not publish"
  )
})

test_that("plan_fetch routes raw and log2 to the requested representation", {
  fx <- make_fixture_cache()

  raw <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    gene_data = fx$genes,
    genes = "G1",
    data_type = "raw",
    offline = TRUE,
    release = fx$version
  )

  log2 <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    gene_data = fx$genes,
    genes = "G1",
    data_type = "log2",
    offline = TRUE,
    release = fx$version
  )

  expect_equal(raw$data_type, "raw")
  expect_true(all(raw$files$data_type == "raw"))
  expect_match(raw$files$relative_path, "-raw\\.h5ad$")

  expect_equal(log2$data_type, "log2")
  expect_true(all(log2$files$data_type == "log2"))
  expect_match(log2$files$relative_path, "-log2\\.h5ad$")
})

test_that("expression matrix reads log2 representation without transforming it", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells
  rownames(cells) <- cells$cell_label

  m <- suppressMessages(
    .abc_expression_matrix(
      abc_cache = cache,
      cells = cells,
      gene_data = fx$genes,
      gene_sel = seq_len(nrow(fx$genes)),
      data_type = "log2"
    )
  )

  expect_s4_class(m, "dgCMatrix")
  expect_equal(
    unname(dense(m)),
    unname(fx$truth_log2),
    tolerance = 1e-12
  )
})

test_that("fetch_data returns raw expression in the default matrix bundle", {

  fx <- make_fixture_cache()

  obj <- suppressMessages(
    fetch_data(
      download_base = fx$dir,
      metadata = fx$cells,
      gene_data = fx$genes,
      genes = c("G1", "G3", "G4"),
      data_type = "raw",
      offline = TRUE,
      release = fx$version
    )
  )

  expect_s3_class(
    obj,
    "abc_matrix"
  )

  expect_identical(
    obj$data_type,
    "raw"
  )

  expect_equal(
    unname(dense(obj$expression)),
    unname(
      fx$truth[
        c("G1", "G3", "G4"),
        ,
        drop = FALSE
      ]
    )
  )
})

test_that("fetch_data returns published log2 expression without transforming it", {

  fx <- make_fixture_cache()

  obj <- suppressMessages(
    fetch_data(
      download_base = fx$dir,
      metadata = fx$cells,
      gene_data = fx$genes,
      genes = c("G1", "G3", "G4"),
      data_type = "log2",
      offline = TRUE,
      release = fx$version
    )
  )

  expect_s3_class(
    obj,
    "abc_matrix"
  )

  expect_identical(
    obj$data_type,
    "log2"
  )

  expect_equal(
    unname(dense(obj$expression)),
    unname(
      fx$truth_log2[
        c("G1", "G3", "G4"),
        ,
        drop = FALSE
      ]
    ),
    tolerance = 1e-12
  )
})

test_that("fetch_data uses the route default representation", {

  fx <- make_fixture_cache()

  obj <- suppressMessages(
    fetch_data(
      download_base = fx$dir,
      metadata = fx$cells,
      gene_data = fx$genes,
      genes = "Ddd",
      offline = TRUE,
      release = fx$version
    )
  )

  expect_s3_class(
    obj,
    "abc_matrix"
  )

  expect_identical(
    obj$data_type,
    "raw"
  )

  expect_identical(
    rownames(obj$expression),
    "G4"
  )
})

test_that("expression inventory reports published representations", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  inventory <- .abc_expression_inventory(cache, "TEST-10X")

  expect_equal(nrow(inventory), 4L)
  expect_setequal(
    paste(
      inventory$feature_matrix_label,
      inventory$data_type,
      sep = "/"
    ),
    c(
      "TEST-A/raw",
      "TEST-A/log2",
      "TEST-B/raw",
      "TEST-B/log2"
    )
  )
})


test_that("route provenance overrides legacy dataset_label routing", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells
  cells$dataset_label <- "NOT-AN-EXPRESSION-DIRECTORY"

  files <- .abc_expression_files(
    cache,
    cells,
    data_type = "raw"
  )

  expect_equal(
    files$dataset_label,
    c("TEST-10X", "TEST-10X")
  )
  expect_equal(
    files$feature_matrix_label,
    c("TEST-A", "TEST-B")
  )
  expect_equal(files$n_cells, c(3L, 3L))
})


test_that("route-aware expression routing rejects unknown matrix labels", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells
  cells$feature_matrix_label[1] <- "OUTSIDE-ROUTE"

  expect_error(
    .abc_expression_files(
      cache,
      cells,
      data_type = "raw"
    ),
    "Cannot route feature_matrix_label='OUTSIDE-ROUTE'"
  )
})


test_that("route-aware expression routing rejects missing matrix labels", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells
  cells$feature_matrix_label[1] <- NA_character_

  expect_error(
    .abc_expression_files(
      cache,
      cells,
      data_type = "raw"
    ),
    "missing feature_matrix_label"
  )
})


test_that("legacy metadata still routes through dataset_label", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells
  attr(cells, "abc_expression_directories") <- NULL

  files <- .abc_expression_files(
    cache,
    cells,
    data_type = "raw"
  )

  expect_equal(
    files$dataset_label,
    c("TEST-10X", "TEST-10X")
  )
  expect_equal(
    files$feature_matrix_label,
    c("TEST-A", "TEST-B")
  )
})

test_that("route-aware expression routing rejects ambiguous matrices", {
  fx <- make_fixture_cache()

  manifest_path <- file.path(
    fx$dir,
    "releases",
    fx$version,
    "manifest.json"
  )

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  manifest$file_listing[["TEST-ALT"]] <-
    manifest$file_listing[["TEST-10X"]]

  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE
  )

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells
  attr(cells, "abc_expression_directories") <-
    c("TEST-10X", "TEST-ALT")

  expect_error(
    .abc_expression_files(
      cache,
      cells,
      data_type = "raw"
    ),
    "Ambiguous expression routing"
  )
})

test_that("expression files outside declared route are not accessible", {
  fx <- make_fixture_cache()

  manifest_path <- file.path(
    fx$dir,
    "releases",
    fx$version,
    "manifest.json"
  )

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  outside <- manifest$file_listing[["TEST-10X"]]
  names(outside$expression_matrices) <-
    paste0(
      names(outside$expression_matrices),
      "-OUTSIDE"
    )

  manifest$file_listing[["TEST-OUTSIDE"]] <- outside

  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE
  )

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells[1, , drop = FALSE]
  cells$feature_matrix_label <- "TEST-A-OUTSIDE"

  # TEST-OUTSIDE really contains this matrix, but this route only permits
  # TEST-10X. Expression discovery must not escape the declared route.
  expect_error(
    .abc_expression_files(
      cache,
      cells,
      data_type = "raw"
    ),
    "Cannot route feature_matrix_label='TEST-A-OUTSIDE'"
  )
})

test_that("gene resolution accepts identifiers and unique symbols automatically", {
  fx <- make_fixture_cache()

  by_id <- .abc_resolve_genes(
    fx$genes,
    c("G4", "G2")
  )

  expect_identical(
    as.character(by_id$gene_identifier),
    c("G4", "G2")
  )
  expect_identical(
    by_id$.abc_gene_index,
    c(4L, 2L)
  )

  by_symbol <- .abc_resolve_genes(
    fx$genes,
    c("Ddd", "Bbb")
  )

  expect_identical(
    as.character(by_symbol$gene_identifier),
    c("G4", "G2")
  )
  expect_identical(
    by_symbol$.abc_gene_index,
    c(4L, 2L)
  )
})

test_that("gene resolution explains ambiguous and missing genes", {
  fx <- make_fixture_cache()

  expect_error(
    .abc_resolve_genes(fx$genes, "Aaa"),
    "G1, G3"
  )

  expect_error(
    .abc_resolve_genes(fx$genes, "NOT-A-GENE"),
    "Gene 'NOT-A-GENE' was not found",
    fixed = TRUE
  )
})

test_that("gene resolution rejects duplicate requests", {
  fx <- make_fixture_cache()

  expect_error(
    .abc_resolve_genes(
      fx$genes,
      c("G2", "G2")
    ),
    "without duplicate requests"
  )
})

test_that("different gene names cannot resolve to the same feature twice", {
  fx <- make_fixture_cache()

  expect_error(
    .abc_resolve_genes(
      fx$genes,
      c("G2", "Bbb")
    ),
    "same gene identifier"
  )
})

test_that("gene resolution supports source identifiers when available", {
  fx <- make_fixture_cache()

  genes <- fx$genes
  genes$source_gene_identifier <- c(
    "SRC1",
    "SRC2",
    "SRC3",
    "SRC4"
  )

  out <- .abc_resolve_genes(
    genes,
    c("SRC3", "SRC1")
  )

  expect_identical(
    as.character(out$gene_identifier),
    c("G3", "G1")
  )
  expect_identical(
    out$.abc_gene_index,
    c(3L, 1L)
  )
})

test_that("fetch_data accepts gene identifiers without extra user options", {

  fx <- make_fixture_cache()

  obj <- suppressMessages(
    fetch_data(
      download_base = fx$dir,
      metadata = fx$cells[1:3, , drop = FALSE],
      gene_data = fx$genes,
      genes = c("G4", "G2"),
      offline = TRUE,
      release = fx$version
    )
  )

  expect_s3_class(
    obj,
    "abc_matrix"
  )

  expect_identical(
    rownames(obj$expression),
    c("G4", "G2")
  )

  expect_identical(
    colnames(obj$expression),
    c("c1", "c2", "c3")
  )

})
test_that("plan_fetch counts resolved identifiers like symbols", {
  fx <- make_fixture_cache()

  plan <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells[1:3, , drop = FALSE],
    gene_data = fx$genes,
    genes = c("G4", "G2"),
    offline = TRUE,
    release = fx$version
  )

  expect_identical(plan$n_genes, 2L)
  expect_identical(plan$n_cells, 3L)
})

test_that("duplicate source genes align to H5AD occurrences one-to-one", {
  gene <- data.frame(
    gene_identifier = c(
      "G1.iabdup1",
      "G2",
      "G1.iabdup2"
    ),
    source_gene_identifier = c("G1", "G2", "G1"),
    source_gene_occurrence = c(1L, 1L, 2L),
    source_gene_group_size = c(2L, 1L, 2L),
    stringsAsFactors = FALSE
  )

  expect_identical(
    .abc_gene_positions(
      gene,
      c(1L, 3L, 2L),
      c("G1", "G2", "G1"),
      "duplicate.h5ad"
    ),
    c(1L, 3L, 2L)
  )
})

test_that("duplicate source alignment rejects mismatched H5AD group sizes", {
  gene <- data.frame(
    gene_identifier = c(
      "G1.iabdup1",
      "G1.iabdup2"
    ),
    source_gene_identifier = c("G1", "G1"),
    source_gene_occurrence = c(1L, 2L),
    source_gene_group_size = c(2L, 2L),
    stringsAsFactors = FALSE
  )

  expect_error(
    .abc_gene_positions(
      gene,
      c(1L, 2L),
      c("G1", "G1", "G1"),
      "bad-duplicate.h5ad"
    ),
    "selected genes are missing"
  )
})

test_that("canonical H5AD identifiers remain preferred", {
  gene <- data.frame(
    gene_identifier = c(
      "G1.iabdup1",
      "G1.iabdup2"
    ),
    source_gene_identifier = c("G1", "G1"),
    source_gene_occurrence = c(1L, 2L),
    source_gene_group_size = c(2L, 2L),
    stringsAsFactors = FALSE
  )

  expect_identical(
    .abc_gene_positions(
      gene,
      c(2L, 1L),
      c("G1.iabdup1", "G1.iabdup2"),
      "canonical.h5ad"
    ),
    c(2L, 1L)
  )
})

test_that("expression matrix aligns normalized duplicate genes to H5AD occurrences", {
  fx <- make_fixture_cache(
    var_ids = c("G1", "G2", "G1", "G4")
  )

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  # Model the real Macaque condition: metadata and H5AD both publish a
  # duplicated stable identifier before metadata normalization.
  source_genes <- fx$genes
  source_genes$gene_identifier <- c("G1", "G2", "G1", "G4")

  genes <- .abc_normalize_gene_metadata(source_genes)

  expect_identical(
    genes$gene_identifier,
    c("G1.iabdup1", "G2", "G1.iabdup2", "G4")
  )

  expect_identical(
    genes$source_gene_identifier,
    c("G1", "G2", "G1", "G4")
  )

  expect_identical(
    genes$source_gene_occurrence,
    c(1L, 1L, 2L, 1L)
  )

  expect_identical(
    genes$source_gene_group_size,
    c(2L, 1L, 2L, 1L)
  )

  cells <- fx$cells
  rownames(cells) <- cells$cell_label

  m <- suppressMessages(
    .abc_expression_matrix(
      cache,
      cells,
      genes,
      c(3L, 1L, 4L)
    )
  )

  expect_s4_class(m, "dgCMatrix")

  # Metadata rows 3, 1 and 4 must resolve to H5AD occurrences 3, 1 and 4.
  expect_equal(
    unname(dense(m)),
    unname(
      fx$truth[
        c("G3", "G1", "G4"),
        ,
        drop = FALSE
      ]
    )
  )

  expect_equal(rownames(m), c("G1.iabdup2", "G1.iabdup1", "G4"))

  expect_equal(
    colnames(m),
    fx$cells$cell_label
  )
})

test_that("expression routing can map published cell labels to another matrix", {
  fx <- make_fixture_cache()

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells[1:3, , drop = FALSE]

  published <- as.character(cells$feature_matrix_label)

  # Route TEST-A cells internally to TEST-B while preserving the published
  # feature_matrix_label column itself.
  attr(cells, "abc_expression_directories") <- "TEST-10X"
  attr(cells, "abc_expression_matrix_map") <- c(
    "TEST-A" = "TEST-B"
  )

  files <- .abc_expression_files(
    cache,
    cells,
    data_type = "raw"
  )

  expect_identical(
    as.character(cells$feature_matrix_label),
    published
  )

  expect_false("TEST-A" %in% files$feature_matrix_label)
  expect_true("TEST-B" %in% files$feature_matrix_label)

  expect_equal(
    sum(files$n_cells),
    nrow(cells)
  )
})

test_that("fetch inputs resolve abc_data provenance", {
  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    stringsAsFactors = FALSE
  )

  gene_data <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  data <- structure(
    list(
      cell_metadata = metadata,
      gene_data = gene_data,
      unique_values = list(),
      dataset = "TEST",
      download_base = "/tmp/test-cache",
      release = "20260711"
    ),
    class = c("abc_data", "list")
  )

  out <- .abc_fetch_inputs(data = data)

  expect_identical(out$metadata, metadata)
  expect_identical(out$gene_data, gene_data)
  expect_identical(
    out$download_base,
    "/tmp/test-cache"
  )
  expect_identical(
    out$release,
    "20260711"
  )
})

test_that("fetch inputs reject conflicting interfaces", {
  data <- structure(
    list(
      cell_metadata = data.frame(cell_label = "c1"),
      gene_data = data.frame(
        gene_identifier = "G1",
        gene_symbol = "A"
      ),
      download_base = "/tmp/test-cache",
      release = "20260711"
    ),
    class = c("abc_data", "list")
  )

  expect_error(
    .abc_fetch_inputs(
      data = data,
      metadata = data$cell_metadata
    ),
    "either data or metadata/gene_data",
    fixed = TRUE
  )
})

test_that("plan_fetch accepts abc_data directly", {
  fx <- make_fixture_cache()

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  built <- suppressMessages(
    .build_dataset(
      fixture_spec,
      cache,
      "TEST"
    )
  )

  data <- structure(
    c(
      built,
      list(
        dataset = "TEST",
        download_base = fx$dir,
        release = fx$version
      )
    ),
    class = c("abc_data", "list")
  )

  plan <- plan_fetch(
    data,
    genes = "G1",
    offline = TRUE
  )

  expect_s3_class(plan, "abc_fetch_plan")
  expect_identical(plan$n_genes, 1L)
  expect_identical(plan$n_cells, 6L)
  expect_identical(plan$release, fx$version)
})

test_that("plan_fetch abc_data and legacy inputs agree", {
  fx <- make_fixture_cache()

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  built <- suppressMessages(
    .build_dataset(
      fixture_spec,
      cache,
      "TEST"
    )
  )

  data <- structure(
    c(
      built,
      list(
        dataset = "TEST",
        download_base = fx$dir,
        release = fx$version
      )
    ),
    class = c("abc_data", "list")
  )

  modern <- plan_fetch(
    data,
    genes = "G1",
    offline = TRUE
  )

  legacy <- plan_fetch(
    metadata = data$cell_metadata,
    gene_data = data$gene_data,
    genes = "G1",
    download_base = fx$dir,
    release = fx$version,
    offline = TRUE
  )

  expect_identical(modern$n_cells, legacy$n_cells)
  expect_identical(modern$n_genes, legacy$n_genes)
  expect_identical(modern$data_type, legacy$data_type)

  expect_identical(
    modern$files$feature_matrix_label,
    legacy$files$feature_matrix_label
  )
})

test_that("load_data output passes directly to plan_fetch", {
  fx <- make_fixture_cache()

  local_mocked_bindings(
    .abc_datasets = list(TEST = fixture_spec)
  )

  data <- suppressMessages(
    load_data(
      "TEST",
      download_base = fx$dir,
      offline = TRUE,
      release = fx$version,
      cell_ids = c("c6", "c1", "c3")
    )
  )

  expect_s3_class(data, "abc_data")
  expect_identical(
    data$cell_metadata$cell_label,
    c("c6", "c1", "c3")
  )

  plan <- plan_fetch(
    data,
    genes = "G1",
    offline = TRUE
  )

  expect_s3_class(plan, "abc_fetch_plan")
  expect_identical(plan$n_cells, 3L)
  expect_identical(plan$n_genes, 1L)
  expect_identical(plan$release, fx$version)
})

test_that("matrix materialization preserves expression and metadata", {
  expression <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 3, 0, 4),
      nrow = 2L,
      byrow = TRUE,
      dimnames = list(c("G1", "G2"), c("c1", "c2", "c3"))
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2", "c3"),
    group = c("X", "Y", "X"),
    stringsAsFactors = FALSE
  )

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "raw",
    return_type = "matrix"
  )

  expect_s3_class(out, "abc_matrix")
  expect_s4_class(out$expression, "dgCMatrix")
  expect_identical(
    names(out$cell_metadata),
    names(metadata)
  )
  expect_identical(
    names(out$gene_data),
    names(genes)
  )

  for (column in names(metadata)) {
    expect_identical(
      out$cell_metadata[[column]],
      metadata[[column]]
    )
  }

  for (column in names(genes)) {
    expect_identical(
      out$gene_data[[column]],
      genes[[column]]
    )
  }

  expect_identical(out$data_type, "raw")

  expect_identical(
    rownames(out$cell_metadata),
    colnames(out$expression)
  )
  expect_identical(
    rownames(out$gene_data),
    rownames(out$expression)
  )
})

test_that("Monocle materialization keeps raw expression and metadata", {

  expression <- Matrix::Matrix(
    matrix(
      c(
        1, 0, 2,
        0, 3, 1
      ),
      nrow = 2L,
      byrow = TRUE,
      dimnames = list(
        c("G1", "G2"),
        c("c1", "c2", "c3")
      )
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2", "c3"),
    cell_type = c("A", "B", "A"),
    stringsAsFactors = FALSE,
    row.names = c("c1", "c2", "c3")
  )

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("GENE_A", "GENE_B"),
    stringsAsFactors = FALSE,
    row.names = c("G1", "G2")
  )

  out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "raw",
    return_type = "monocle",
    provenance = list(
      dataset = "TEST"
    )
  )

  expect_s4_class(
    out,
    "CellDataSet"
  )

  expect_equal(
    unname(as.matrix(Biobase::exprs(out))),
    unname(as.matrix(expression))
  )

  expect_identical(
    rownames(Biobase::pData(out)),
    c("c1", "c2", "c3")
  )

  expect_identical(
    Biobase::pData(out)$cell_type,
    c("A", "B", "A")
  )

  expect_identical(
    rownames(Biobase::fData(out)),
    c("G1", "G2")
  )

  expect_identical(
    Biobase::fData(out)$gene_short_name,
    c("GENE_A", "GENE_B")
  )

  expect_identical(
    Biobase::experimentData(out)@other$ABCatlas$dataset,
    "TEST"
  )
})

test_that("Monocle materialization rejects published log2 expression", {

  expression <- Matrix::Matrix(
    matrix(
      c(
        0.5, 0,
        1.5, 2
      ),
      nrow = 2L,
      byrow = TRUE,
      dimnames = list(
        c("G1", "G2"),
        c("c1", "c2")
      )
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    stringsAsFactors = FALSE,
    row.names = c("c1", "c2")
  )

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("GENE_A", "GENE_B"),
    stringsAsFactors = FALSE,
    row.names = c("G1", "G2")
  )

  expect_error(
    .abc_materialize(
      expression,
      metadata,
      genes,
      data_type = "log2",
      return_type = "monocle"
    ),
    "requires raw expression data"
  )
})

test_that("spatial coordinate selection prefers canonical route coordinates", {
  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    x_experiment = c(100, 200),
    y_experiment = c(300, 400),
    x_rotated_um = c(10, 20),
    y_rotated_um = c(30, 40),
    x_slab_mm = c(1, 2),
    y_slab_mm = c(3, 4),
    stringsAsFactors = FALSE
  )

  coords <- .abc_spatial_coordinates(metadata)

  expect_identical(
    colnames(coords),
    c("x_slab_mm", "y_slab_mm")
  )

  expect_equal(
    unname(coords),
    matrix(c(1, 2, 3, 4), ncol = 2L)
  )
})

test_that("spatial coordinate selection accepts an explicit pair", {
  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    x_experiment = c(100, 200),
    y_experiment = c(300, 400),
    x_slab_mm = c(1, 2),
    y_slab_mm = c(3, 4),
    stringsAsFactors = FALSE
  )

  coords <- .abc_spatial_coordinates(
    metadata,
    coordinate_columns = c(
      "x_experiment",
      "y_experiment"
    )
  )

  expect_identical(
    colnames(coords),
    c("x_experiment", "y_experiment")
  )
})

test_that("SingleCellExperiment raw materialization preserves count semantics", {
  expression <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 3, 0, 4),
      nrow = 2L,
      dimnames = list(c("G1", "G2"), c("c1", "c2", "c3"))
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2", "c3"),
    group = c("X", "Y", "X"),
    stringsAsFactors = FALSE
  )

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "raw",
    return_type = "sce"
  )

  expect_s4_class(out, "SingleCellExperiment")
  expect_identical(
    SummarizedExperiment::assayNames(out),
    "counts"
  )

  expect_equal(
    unname(as.matrix(
      SummarizedExperiment::assay(out, "counts")
    )),
    unname(as.matrix(expression))
  )

  expect_identical(
    colnames(out),
    metadata$cell_label
  )
  expect_identical(
    rownames(out),
    rownames(expression)
  )

  expect_identical(
    as.character(
      SummarizedExperiment::colData(out)$cell_label
    ),
    metadata$cell_label
  )

  expect_identical(
    as.character(
      SummarizedExperiment::rowData(out)$gene_identifier
    ),
    genes$gene_identifier
  )
})

test_that("SingleCellExperiment log2 materialization does not invent counts", {
  expression <- Matrix::Matrix(
    matrix(
      c(0.5, 0, 1.5, 2),
      nrow = 2L,
      dimnames = list(c("G1", "G2"), c("c1", "c2"))
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    stringsAsFactors = FALSE
  )

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "log2",
    return_type = "sce"
  )

  expect_s4_class(out, "SingleCellExperiment")

  expect_identical(
    SummarizedExperiment::assayNames(out),
    "log2"
  )

  expect_false(
    "counts" %in%
      SummarizedExperiment::assayNames(out)
  )

  expect_equal(
    unname(as.matrix(
      SummarizedExperiment::assay(out, "log2")
    )),
    unname(as.matrix(expression))
  )
})

test_that("SpatialExperiment raw materialization preserves coordinates and counts", {
  expression <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 3),
      nrow = 2L,
      dimnames = list(c("G1", "G2"), c("c1", "c2"))
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    x_experiment = c(100, 200),
    y_experiment = c(300, 400),
    x_slab_mm = c(1, 2),
    y_slab_mm = c(3, 4),
    stringsAsFactors = FALSE
  )

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "raw",
    return_type = "spatial",
    spatial = TRUE
  )

  expect_s4_class(out, "SpatialExperiment")

  expect_identical(
    SummarizedExperiment::assayNames(out),
    "counts"
  )

  expect_equal(
    unname(as.matrix(
      SummarizedExperiment::assay(out, "counts")
    )),
    unname(as.matrix(expression))
  )

  coords <- SpatialExperiment::spatialCoords(out)

  expect_identical(
    colnames(coords),
    c("x_slab_mm", "y_slab_mm")
  )

  expect_equal(
    unname(coords),
    matrix(
      c(1, 2, 3, 4),
      ncol = 2L
    )
  )

  expect_identical(
    rownames(coords),
    metadata$cell_label
  )
})

test_that("SpatialExperiment log2 materialization preserves transformed semantics", {
  expression <- Matrix::Matrix(
    matrix(
      c(0.5, 0, 1.5, 2),
      nrow = 2L,
      dimnames = list(c("G1", "G2"), c("c1", "c2"))
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    x_experiment = c(100, 200),
    y_experiment = c(300, 400),
    stringsAsFactors = FALSE
  )

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "log2",
    return_type = "spatial",
    coordinate_columns = c(
      "x_experiment",
      "y_experiment"
    )
  )

  expect_s4_class(out, "SpatialExperiment")

  expect_identical(
    SummarizedExperiment::assayNames(out),
    "log2"
  )

  expect_false(
    "counts" %in%
      SummarizedExperiment::assayNames(out)
  )

  expect_equal(
    unname(as.matrix(
      SummarizedExperiment::assay(out, "log2")
    )),
    unname(as.matrix(expression))
  )

  expect_identical(
    colnames(
      SpatialExperiment::spatialCoords(out)
    ),
    c("x_experiment", "y_experiment")
  )
})

test_that("spatial output requires declared modality or explicit coordinates", {
  expression <- Matrix::Matrix(
    matrix(
      c(1, 2),
      nrow = 1L,
      dimnames = list(
        "G1",
        c("c1", "c2")
      )
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    x_experiment = c(100, 200),
    y_experiment = c(300, 400),
    stringsAsFactors = FALSE
  )

  genes <- data.frame(
    gene_identifier = "G1",
    gene_symbol = "A",
    stringsAsFactors = FALSE
  )

  expect_error(
    .abc_materialize(
      expression,
      metadata,
      genes,
      data_type = "raw",
      return_type = "spatial"
    ),
    "No declared spatial modality",
    fixed = TRUE
  )

  out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "raw",
    return_type = "spatial",
    coordinate_columns = c(
      "x_experiment",
      "y_experiment"
    )
  )

  expect_s4_class(out, "SpatialExperiment")

  expect_identical(
    colnames(
      SpatialExperiment::spatialCoords(out)
    ),
    c("x_experiment", "y_experiment")
  )
})

test_that("all return formats preserve ABC provenance", {
  expression <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 3),
      nrow = 2L,
      dimnames = list(
        c("G1", "G2"),
        c("c1", "c2")
      )
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    x_experiment = c(100, 200),
    y_experiment = c(300, 400),
    stringsAsFactors = FALSE
  )

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  provenance <- list(
    dataset = "TEST",
    release = "20260711",
    data_type = "raw",
    modality = "spatial_transcriptomics",
    spatial = TRUE,
    expression_directories = "TEST-DIR",
    n_cells = 2L,
    n_genes = 2L
  )

  matrix_out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "raw",
    return_type = "matrix",
    provenance = provenance,
    spatial = TRUE
  )

  expect_identical(
    matrix_out$provenance,
    provenance
  )

  monocle_out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "raw",
    return_type = "monocle",
    provenance = provenance,
    spatial = TRUE
  )

  expect_s4_class(
    monocle_out,
    "CellDataSet"
  )

  expect_identical(
    Biobase::experimentData(
      monocle_out
    )@other$ABCatlas,
    provenance
  )

  sce_out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "raw",
    return_type = "sce",
    provenance = provenance,
    spatial = TRUE
  )

  expect_identical(
    S4Vectors::metadata(sce_out)$ABCatlas,
    provenance
  )

  spatial_out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "raw",
    return_type = "spatial",
    provenance = provenance,
    spatial = TRUE
  )

  spatial_provenance <-
    S4Vectors::metadata(spatial_out)$ABCatlas

  expect_identical(
    spatial_provenance$dataset,
    "TEST"
  )
  expect_identical(
    spatial_provenance$release,
    "20260711"
  )
  expect_identical(
    spatial_provenance$data_type,
    "raw"
  )
  expect_identical(
    spatial_provenance$modality,
    "spatial_transcriptomics"
  )
  expect_true(spatial_provenance$spatial)

  expect_identical(
    spatial_provenance$coordinate_columns,
    c("x_experiment", "y_experiment")
  )
})

test_that("native output provenance records route and selection context", {
  metadata <- data.frame(
    cell_label = c("c1", "c2", "c3"),
    stringsAsFactors = FALSE
  )

  attr(metadata, "abc_dataset") <- "TEST"
  attr(metadata, "abc_modality") <- "imputed_spatial_transcriptomics"
  attr(metadata, "abc_spatial") <- TRUE
  attr(metadata, "abc_expression_directories") <- "TEST-IMPUTED"

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  provenance <- .abc_output_provenance(
    metadata,
    genes,
    data_type = "log2",
    release = "20260711"
  )

  expect_identical(provenance$dataset, "TEST")
  expect_identical(provenance$release, "20260711")
  expect_identical(provenance$data_type, "log2")
  expect_identical(
    provenance$modality,
    "imputed_spatial_transcriptomics"
  )
  expect_true(provenance$spatial)

  expect_identical(
    provenance$expression_directories,
    "TEST-IMPUTED"
  )

  expect_identical(provenance$n_cells, 3L)
  expect_identical(provenance$n_genes, 2L)
})

test_that("matrix materialization does not require optional output packages", {
  expression <- Matrix::Matrix(
    matrix(
      c(1, 2),
      nrow = 1L,
      dimnames = list(
        "G1",
        c("c1", "c2")
      )
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    stringsAsFactors = FALSE
  )

  genes <- data.frame(
    gene_identifier = "G1",
    gene_symbol = "A",
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    .abc_require_output_package = function(...) {
      stop(
        "optional package lookup should not occur",
        call. = FALSE
      )
    }
  )

  out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "raw",
    return_type = "matrix"
  )

  expect_s3_class(out, "abc_matrix")
  expect_s4_class(out$expression, "dgCMatrix")
})

test_that("H5AD expression reader supports sparse and dense X layouts", {
  fx <- make_fixture_cache()

  sparse_path <- file.path(
    fx$dir,
    "expression_matrices",
    "TEST-10X",
    fx$version,
    "TEST-A-raw.h5ad"
  )

  sparse <- .abc_h5ad_expression(
    sparse_path,
    n_genes = 4L,
    n_cells = 3L
  )

  expect_s4_class(
    sparse,
    "H5SparseMatrix"
  )

  dense_path <- file.path(
    withr::local_tempdir(),
    "dense.h5ad"
  )

  counts <- matrix(
    c(
      1, 0, 5, 0,
      0, 2, 0, 7,
      3, 0, 0, 4
    ),
    nrow = 3L,
    byrow = TRUE
  )

  write_dense_h5ad(
    dense_path,
    counts = counts,
    cell_labels = c("c1", "c2", "c3"),
    gene_ids = c("G1", "G2", "G3", "G4"),
    gene_symbols = c("A", "B", "C", "D")
  )

  dense_x <- .abc_h5ad_expression(
    dense_path,
    n_genes = 4L,
    n_cells = 3L
  )

  expect_s4_class(
    dense_x,
    "HDF5Matrix"
  )

  expect_identical(
    as.integer(dim(dense_x)),
    c(4L, 3L)
  )

  expect_equal(
    unname(as.matrix(dense_x)),
    unname(t(counts))
  )
})


test_that("dense H5AD expression preserves requested gene and cell order", {
  dir <- withr::local_tempdir()
  path <- file.path(dir, "dense.h5ad")

  cells <- c("c1", "c2", "c3")

  genes <- data.frame(
    gene_identifier = c(
      "G1",
      "G2",
      "G3",
      "G4"
    ),
    gene_symbol = c(
      "A",
      "B",
      "C",
      "D"
    ),
    stringsAsFactors = FALSE
  )

  # cells x genes
  counts <- matrix(
    c(
      1, 0, 5, 0,
      0, 2, 0, 7,
      3, 0, 0, 4
    ),
    nrow = 3L,
    byrow = TRUE
  )

  write_dense_h5ad(
    path,
    counts = counts,
    cell_labels = cells,
    gene_ids = genes$gene_identifier,
    gene_symbols = genes$gene_symbol
  )

  var <- .h5ad_index(
    path,
    "var"
  )
  obs <- .h5ad_index(
    path,
    "obs"
  )

  gene_sel <- c(4L, 1L, 3L)
  requested_cells <- c(
    "c3",
    "c1"
  )

  gene_pos <- .abc_gene_positions(
    genes,
    gene_sel,
    var,
    basename(path)
  )

  cell_pos <- match(
    requested_cells,
    obs
  )

  X <- .abc_h5ad_expression(
    path,
    n_genes = length(var),
    n_cells = length(obs)
  )

  out <- methods::as(
    X[
      gene_pos,
      cell_pos,
      drop = FALSE
    ],
    "dgCMatrix"
  )

  expect_identical(
    dim(out),
    c(3L, 2L)
  )

  expect_equal(
    unname(as.matrix(out)),
    unname(
      t(counts)[
        c(4L, 1L, 3L),
        c(3L, 1L),
        drop = FALSE
      ]
    )
  )
})


test_that("H5AD expression reader rejects inconsistent dimensions", {
  path <- file.path(
    withr::local_tempdir(),
    "dense-bad-dimensions.h5ad"
  )

  counts <- matrix(
    seq_len(12L),
    nrow = 3L,
    ncol = 4L
  )

  write_dense_h5ad(
    path,
    counts = counts,
    cell_labels = c("c1", "c2", "c3"),
    gene_ids = c("G1", "G2", "G3", "G4"),
    gene_symbols = c("A", "B", "C", "D")
  )

  expect_error(
    .abc_h5ad_expression(
      path,
      n_genes = 3L,
      n_cells = 4L
    ),
    "/X dimensions are",
    fixed = TRUE
  )

  expect_error(
    .abc_h5ad_expression(
      path,
      n_genes = 3L,
      n_cells = 4L
    ),
    "Automatic transposition is intentionally disabled",
    fixed = TRUE
  )
})

test_that("accession-qualified NCBI H5AD identifiers normalize safely", {
  input <- c(
    "NCBIGene:114678393:NC_041754.1",
    "NCBIGene:114675313:NW_021160097.1",
    "NCBIGene:710010",
    "ENSMUSG00000000001",
    "NCBIGene:12345:OTHER_123.1"
  )

  expect_identical(
    .abc_h5ad_source_gene_identifier(input),
    c(
      "NCBIGene:114678393",
      "NCBIGene:114675313",
      "NCBIGene:710010",
      "ENSMUSG00000000001",
      "NCBIGene:12345:OTHER_123.1"
    )
  )
})


test_that("gene positions resolve accession-qualified NCBI features", {
  gene_data <- data.frame(
    gene_identifier = c(
      "NCBIGene:114678393",
      "NCBIGene:114675313"
    ),
    source_gene_identifier = c(
      "NCBIGene:114678393",
      "NCBIGene:114675313"
    ),
    source_gene_occurrence = c(1L, 1L),
    source_gene_group_size = c(1L, 1L),
    stringsAsFactors = FALSE
  )

  var <- c(
    "NCBIGene:114675313:NW_021160097.1",
    "NCBIGene:114678393:NC_041754.1"
  )

  expect_identical(
    .abc_gene_positions(
      gene_data,
      sel = c(1L, 2L),
      var = var,
      file = "TEST-MACAQUE.h5ad"
    ),
    c(2L, 1L)
  )
})


test_that("accession-qualified duplicate genes retain occurrence identity", {
  gene_data <- data.frame(
    gene_identifier = c(
      "NCBIGene:693935.iabdup1",
      "NCBIGene:693935.iabdup2"
    ),
    source_gene_identifier = c(
      "NCBIGene:693935",
      "NCBIGene:693935"
    ),
    source_gene_occurrence = c(1L, 2L),
    source_gene_group_size = c(2L, 2L),
    stringsAsFactors = FALSE
  )

  var <- c(
    "NCBIGene:693935:NC_041754.1",
    "NCBIGene:693935:NC_041754.1"
  )

  expect_identical(
    .abc_gene_positions(
      gene_data,
      sel = c(2L, 1L),
      var = var,
      file = "TEST-MACAQUE-DUPLICATE.h5ad"
    ),
    c(2L, 1L)
  )
})


test_that("unrecognized H5AD identifier suffixes are not stripped", {
  gene_data <- data.frame(
    gene_identifier = "NCBIGene:12345",
    source_gene_identifier = "NCBIGene:12345",
    source_gene_occurrence = 1L,
    source_gene_group_size = 1L,
    stringsAsFactors = FALSE
  )

  var <- c(
    "NCBIGene:12345:OTHER_123.1",
    "NCBIGene:99999"
  )

  expect_identical(
    .abc_h5ad_source_gene_identifier(var),
    var
  )

  expect_error(
    .abc_gene_positions(
      gene_data,
      sel = 1L,
      var = var,
      file = "TEST-UNRECOGNIZED.h5ad"
    ),
    "missing from TEST-UNRECOGNIZED.h5ad",
    fixed = TRUE
  )
})

test_that("materialization validates expression and metadata dimensions", {
  expression <- Matrix::Matrix(
    matrix(
      c(1, 2, 3, 4),
      nrow = 2L,
      dimnames = list(
        c("G1", "G2"),
        c("c1", "c2")
      )
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    stringsAsFactors = FALSE
  )

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    .abc_materialize(
      expression[, 1, drop = FALSE],
      metadata,
      genes,
      data_type = "raw",
      return_type = "matrix"
    ),
    "Expression columns do not match the selected cell metadata.",
    fixed = TRUE
  )

  expect_error(
    .abc_materialize(
      expression[1, , drop = FALSE],
      metadata,
      genes,
      data_type = "raw",
      return_type = "matrix"
    ),
    "Expression rows do not match the selected gene metadata.",
    fixed = TRUE
  )
})

test_that("materialization validates cell alignment", {
  expression <- Matrix::Matrix(
    matrix(
      c(1, 2, 3, 4),
      nrow = 2L,
      dimnames = list(
        c("G1", "G2"),
        c("c2", "c1")
      )
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    stringsAsFactors = FALSE
  )

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    .abc_materialize(
      expression,
      metadata,
      genes,
      data_type = "raw",
      return_type = "matrix"
    ),
    "Expression columns are not aligned to cell_metadata$cell_label.",
    fixed = TRUE
  )
})

test_that("materialization requires unique non-missing gene identifiers", {
  expression <- Matrix::Matrix(
    matrix(
      c(1, 2, 3, 4),
      nrow = 2L,
      dimnames = list(
        c("G1", "G2"),
        c("c1", "c2")
      )
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    stringsAsFactors = FALSE
  )

  genes <- data.frame(
    gene_identifier = c("G1", "G1"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    .abc_materialize(
      expression,
      metadata,
      genes,
      data_type = "raw",
      return_type = "matrix"
    ),
    "gene_data$gene_identifier must be a non-missing unique key.",
    fixed = TRUE
  )
})

test_that("materialization validates gene alignment", {
  expression <- Matrix::Matrix(
    matrix(
      c(1, 2, 3, 4),
      nrow = 2L,
      dimnames = list(
        c("G2", "G1"),
        c("c1", "c2")
      )
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    stringsAsFactors = FALSE
  )

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    .abc_materialize(
      expression,
      metadata,
      genes,
      data_type = "raw",
      return_type = "matrix"
    ),
    "Expression rows are not aligned to gene_data$gene_identifier.",
    fixed = TRUE
  )
})

test_that("materialization builds default provenance from metadata", {
  expression <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 3),
      nrow = 2L,
      dimnames = list(
        c("G1", "G2"),
        c("c1", "c2")
      )
    ),
    sparse = TRUE
  )

  metadata <- data.frame(
    cell_label = c("c1", "c2"),
    stringsAsFactors = FALSE
  )

  attr(metadata, "abc_modality") <- "single_cell_transcriptomics"
  attr(metadata, "abc_spatial") <- FALSE

  genes <- data.frame(
    gene_identifier = c("G1", "G2"),
    gene_symbol = c("A", "B"),
    stringsAsFactors = FALSE
  )

  out <- .abc_materialize(
    expression,
    metadata,
    genes,
    data_type = "log2",
    return_type = "matrix"
  )

  expect_identical(
    out$provenance,
    list(
      data_type = "log2",
      modality = "single_cell_transcriptomics",
      spatial = FALSE,
      n_cells = 2L,
      n_genes = 2L
    )
  )
})

test_that("invalid output requests fail before expression is read", {
  fx <- make_fixture_cache()

  local_mocked_bindings(
    .abc_expression_matrix = function(...) {
      stop("expression matrix should not be read", call. = FALSE)
    }
  )

  base_args <- list(
    download_base = fx$dir,
    metadata = fx$cells,
    gene_data = fx$genes,
    genes = "G1",
    offline = TRUE,
    release = fx$version
  )

  expect_error(
    do.call(fetch_data, c(base_args, list(return_type = "banana"))),
    "should be one of"
  )

  expect_error(
    suppressMessages(do.call(
      fetch_data,
      c(base_args, list(return_type = "monocle", data_type = "log2"))
    )),
    "requires raw expression data"
  )

  expect_error(
    suppressMessages(do.call(
      fetch_data,
      c(base_args, list(return_type = "spatial"))
    )),
    "No declared spatial modality"
  )

  expect_error(
    suppressMessages(do.call(
      fetch_data,
      c(base_args, list(
        return_type = "spatial",
        coordinate_columns = c("x_missing", "y_missing")
      ))
    )),
    "not found in cell metadata"
  )
})

test_that("a plan cannot be combined with arguments it would ignore", {
  fx <- make_fixture_cache()

  plan <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    gene_data = fx$genes,
    genes = "G1",
    offline = TRUE,
    release = fx$version,
    confirm = FALSE
  )

  expect_error(
    fetch_data(plan = plan, genes = "G4"),
    "Remove: genes"
  )

  expect_error(
    fetch_data(
      plan = plan,
      data_type = "log2",
      filters = list(donor_sex = "M")
    ),
    "Remove: filters, data_type"
  )
})

test_that("plan printing aligns the storage lines in both branches", {
  example_dir <- system.file("extdata", "ex", package = "ABCatlas")
  example_data <- readRDS(file.path(example_dir, "example-data.rds"))
  plan <- suppressMessages(plan_fetch(
    download_base = example_dir,
    metadata = example_data$cells,
    gene_data = example_data$genes,
    genes = c("G1", "G4"),
    offline = TRUE,
    release = "20250101",
    confirm = FALSE
  ))

  plan$available_space <- NA_real_
  out <- capture.output(print(plan))
  expect_true("  available disk space: unable to determine" %in% out)
  expect_true("  storage check:        unable to determine" %in% out)

  plan$available_space <- 1e12
  plan$storage_sufficient <- TRUE
  out <- capture.output(print(plan))
  expect_true(any(grepl("^  available disk space: [0-9.]+ [kMGT]?B$", out)))
  expect_true("  storage check:        PASS" %in% out)
})
