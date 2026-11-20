# Behaviour of branches that the main test files do not exercise offline.
# These tests were written against the pre-refactor code and must keep passing
# after the functions involved are split into helpers.

example_dir <- function() {
    system.file("extdata", "ex", package = "ABCatlasR")
}

s3_listing_page <- function(prefixes, truncated, token = NULL) {
    xml <- paste0(
        "<ListBucketResult xmlns=\"http://s3.amazonaws.com/doc/2006-03-01/\">",
        paste0(
            "<CommonPrefixes><Prefix>", prefixes, "</Prefix></CommonPrefixes>",
            collapse = ""
        ),
        "<IsTruncated>", truncated, "</IsTruncated>",
        if (!is.null(token)) {
            paste0("<NextContinuationToken>", token, "</NextContinuationToken>")
        },
        "</ListBucketResult>"
    )
    list(status_code = 200L, content = charToRaw(xml))
}

test_that("remote manifest listing follows S3 continuation pages", {
    urls <- character()
    local_mocked_bindings(
        .abc_fetch_memory_retry = function(url, ...) {
            urls <<- c(urls, url)
            if (!grepl("continuation-token", url, fixed = TRUE)) {
                s3_listing_page(
                    c("releases/20250101/", "releases/20260711/"), "true", "tok/1"
                )
            } else {
                s3_listing_page(
                    c("releases/20240101/", "releases/20250101/"), "false"
                )
            }
        }
    )

    expect_identical(
        .abc_remote_manifests(),
        c(
            "releases/20240101/manifest.json",
            "releases/20250101/manifest.json",
            "releases/20260711/manifest.json"
        )
    )
    expect_length(urls, 2L)
    expect_match(urls[[2L]], "continuation-token=tok%2F1", fixed = TRUE)
})

test_that("remote manifest listing reports HTTP errors", {
    local_mocked_bindings(
        .abc_fetch_memory_retry = function(url, ...) {
            list(status_code = 503L, content = raw())
        }
    )
    expect_error(.abc_remote_manifests(), "HTTP 503", fixed = TRUE)
})

test_that("plan_fetch uses the package release when data records none", {
    abc <- example_abc_data()
    abc$release <- NULL
    expect_error(
        plan_fetch(abc, offline = TRUE),
        paste0("'", .abc_default_release, "' is not available"),
        fixed = TRUE
    )
})

test_that("release = NULL selects the latest local manifest offline", {
    cache <- abc_cache_connect(example_dir(), offline = TRUE, release = NULL)
    expect_identical(cache$manifest_name, "releases/20250101/manifest.json")
    expect_s3_class(cache, "abc_cache")
})

test_that("spatial coordinates reject unusable inputs", {
    md <- data.frame(
        cell_label = c("a", "b"),
        x = c(1, 2),
        y = c(3, NA),
        u = 1:2
    )
    expect_error(
        .abc_spatial_coordinates(md[, c("cell_label", "u")]),
        "No supported spatial coordinate pair",
        fixed = TRUE
    )
    expect_error(
        .abc_spatial_coordinates(md),
        "finite and non-missing",
        fixed = TRUE
    )
})

test_that("catalog rows fall back to defaults for sparse registry entries", {
    row <- .abc_catalog_row("Bare", list())
    expect_identical(row$available_data_types, "raw,log2")
    expect_identical(row$default_data_type, "raw")
    expect_true(is.na(row$cell_directory) && is.na(row$gene_directory))
    expect_true(is.na(row$taxonomy))
    expect_identical(row$expression_directories, "")
    expect_identical(row$modality, "transcriptomics")
})

test_that("files served from the cache are reported", {
    abc <- example_abc_data()
    expect_message(
        fetch_data(
            abc,
            cell_ids = abc$cell_metadata$cell_label,
            genes = c("G1", "G4"),
            offline = TRUE
        ),
        "Using cached .*no download needed"
    )
})
