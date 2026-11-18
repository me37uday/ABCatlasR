# =============================================================================
# ABC Atlas cache: release manifests, file downloads and metadata tables.
#
# Re-implements the parts of abc_atlas_access's AbcProjectCache that the package
# uses, with the same on-disk layout so a download directory can be shared with
# the Python package (paths relative to <download_base>):
#   releases/<version>/manifest.json   release manifests
#   <relative_path>                    data files, as listed in the manifest
#   _manifest_last_used.txt            manifest pinned for this cache
#   _downloaded_data.json              files fully downloaded and verified
# =============================================================================

# Public S3 bucket holding every ABC Atlas release.
.abc_bucket_url <- "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"

# Manifest versions the loaders understand (same window as abc_atlas_access).
.abc_manifest_compatibility <- c("20230101", "20300101")
# ABC Atlas release validated against this package's dataset registry.
.abc_default_release <- "20260711"
# Strings pandas.read_csv treats as missing by default; reused so metadata
# tables have the same missing values as with the Python package.
.abc_na_strings <- c(
    "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN",
    "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN",
    "None", "n/a", "nan", "null"
)

#' Connect to the ABC Atlas project cache
#'
#' Selects a release manifest and returns a handle used to locate, download and
#' read ABC Atlas files under \code{download_base}. With \code{offline = FALSE}
#' (default) the available releases are listed from S3 and missing files are
#' downloaded (and md5-verified) on first use; with \code{offline = TRUE} only
#' manifests and files already on disk are used, and the network is never
#' contacted and nothing is written to \code{download_base}.
#'
#' By default the release validated against this package is used. With
#' \code{release = NULL}, the manifest is chosen as in the Python
#' \code{abc_atlas_access} package: the one last used with this directory if
#' any (with a warning when a newer release exists), otherwise the latest
#' release. The directory layout is the same as in the Python package, so a
#' download directory can be shared between them.
#'
#' @param download_base Download/cache directory; resolved by
#'   \code{\link{abc_download_dir}}.
#' @param offline If \code{TRUE}, use only locally cached data and never contact
#'   S3. Requires a manifest to have been downloaded previously.
#' @param release Optional ABC Atlas release version in \code{YYYYMMDD} form,
#'   for example \code{"20260711"}. When supplied, exactly that release is used.
#'   When \code{NULL}, the previously used/latest manifest selection behavior is
#'   preserved.
#' @return An \code{abc_cache} object: a list with the \code{download_base},
#'   \code{manifest_name}, \code{version} and \code{offline} fields, plus the
#'   functions \code{get_file_path(directory, file_name)} (returns the local
#'   path of a file, downloading it if needed) and
#'   \code{get_metadata_dataframe(directory, file_name, col_classes, ...)}
#'   (reads a metadata CSV into a data.frame; \code{col_classes} defaults to
#'   \code{NULL} and \code{...} goes to \code{data.table::fread}).
#' @examples
#' example_dir <- system.file(
#'     "extdata", "ex",
#'     package = "ABCatlas"
#' )
#'
#' cache <- abc_cache_connect(
#'     example_dir,
#'     offline = TRUE,
#'     release = "20250101"
#' )
#' cache
#'
#' @export
abc_cache_connect <- function(
        download_base = NULL, offline = FALSE,
        release = .abc_default_release) {
    download_base <- abc_download_dir(download_base)
    offline <- isTRUE(offline)

    # Validate offline preconditions up front so an offline call with no local
    # data fails fast with an actionable message.
    if (offline && !.abc_manifest_exists(download_base)) {
        stop(
            "offline = TRUE but no local manifest was found under '",
            download_base,
            "'. Run once with offline = FALSE to download the manifest ",
            "and data."
        )
    }

    available <- if (offline) {
        .abc_local_manifests(download_base)
    } else {
        .abc_remote_manifests()
    }
    if (length(available) == 0L) {
        stop("No ABC Atlas release manifests found.")
    }

    manifest_name <- if (is.null(release)) {
        .abc_choose_manifest(download_base, available)
    } else {
        .abc_release_manifest(release, available, offline)
    }

    manifest_path <- file.path(download_base, manifest_name)
    if (!file.exists(manifest_path)) {
        .abc_download(
            paste0(.abc_bucket_url, "/", manifest_name), manifest_path
        )
    }
    manifest <- jsonlite::read_json(manifest_path)
    .abc_check_manifest_version(manifest$version)

    # Offline connections only read local data and never modify download_base.
    if (!offline) {
        .abc_pin_manifest(download_base, manifest_name)
    }

    .abc_new_cache(download_base, manifest_name, manifest, offline)
}

# Manifest name for an explicitly requested release: a single YYYYMMDD
# version that must be among the available manifests.
.abc_release_manifest <- function(release, available, offline) {
    release <- as.character(release)
    if (length(release) != 1L || is.na(release) ||
        !grepl("^[0-9]{8}$", release)) {
        stop(
            "'release' must be a single date-like version such as '20260711'.",
            call. = FALSE
        )
    }

    manifest_name <- paste0("releases/", release, "/manifest.json")
    if (!manifest_name %in% available) {
        stop(
            "ABC Atlas release '", release, "' is not available",
            if (offline) " in this offline cache" else "", ".",
            call. = FALSE
        )
    }

    manifest_name
}

# Pin the manifest for later sessions. Written without a trailing newline, as
# abc_atlas_access compares the file content verbatim. A failed write is
# reported as a warning and does not stop the connection.
.abc_pin_manifest <- function(download_base, manifest_name) {
    pin_failed <- function(cond) {
        warning(
            "Could not record the manifest in use: ", conditionMessage(cond),
            call. = FALSE
        )
    }

    tryCatch(
        writeLines(
            manifest_name,
            file.path(download_base, "_manifest_last_used.txt"),
            sep = ""
        ),
        warning = pin_failed,
        error = pin_failed
    )
    invisible(NULL)
}

# The abc_cache object: connection fields plus closures that locate,
# download and read files listed in the manifest.
.abc_new_cache <- function(download_base, manifest_name, manifest, offline) {
    cache <- list(
        download_base = download_base,
        manifest_name = manifest_name,
        version = manifest$version,
        offline = offline,
        manifest = manifest
    )
    cache$get_file_path <- function(
            directory, file_name,
            force_download = FALSE,
            skip_hash_check = FALSE) {
        .abc_get_file_path(
            cache, directory, file_name, force_download, skip_hash_check
        )
    }
    cache$get_metadata_dataframe <- function(
            directory, file_name,
            col_classes = NULL, ...) {
        .abc_read_csv(
            cache$get_file_path(directory, file_name), col_classes, ...
        )
    }
    structure(cache, class = "abc_cache")
}

#' @export
print.abc_cache <- function(x, ...) {
    cat("<abc_cache> ", x$manifest_name, if (x$offline) " (offline)", "\n",
        "  download_base: ", x$download_base, "\n",
        "  directories:   ", length(x$manifest$file_listing), "\n",
        sep = ""
    )
    invisible(x)
}

# ---- Manifests --------------------------------------------------------------

# TRUE if a release manifest has already been downloaded into download_base.
.abc_manifest_exists <- function(download_base) {
    length(.abc_local_manifests(download_base)) > 0L
}

# Downloaded manifests as sorted "releases/<version>/manifest.json" names.
.abc_local_manifests <- function(download_base) {
    found <- Sys.glob(
        file.path(download_base, "releases", "*", "manifest.json")
    )
    if (length(found) == 0L) {
        return(character())
    }
    sort(paste0("releases/", basename(dirname(found)), "/manifest.json"),
        method = "radix"
    )
}

.abc_fetch_memory_retry <- function(
        url,
        max_tries = 5L,
        connect_timeout = 10,
        timeout = 30
) {
    max_tries <- as.integer(max_tries)

    if (is.na(max_tries) || max_tries < 1L) {
        stop(
            "`max_tries` must be a positive integer.",
            call. = FALSE
        )
    }

    last_error <- NULL

    for (attempt in seq_len(max_tries)) {
        handle <- curl::new_handle()

        curl::handle_setopt(
            handle,
            connecttimeout = connect_timeout,
            timeout = timeout
        )

        result <- tryCatch(
            curl::curl_fetch_memory(
                url,
                handle = handle
            ),
            error = function(e) {
                e
            }
        )

        if (!inherits(result, "error")) {
            return(result)
        }

        last_error <- result

        if (attempt < max_tries) {
            message(
                "Request attempt ",
                attempt,
                " failed: ",
                conditionMessage(result),
                ". Retrying."
            )
        }
    }

    stop(
        "Request failed after ",
        max_tries,
        " attempts: ",
        conditionMessage(last_error),
        call. = FALSE
    )
}

# Manifests published on S3 (every .json key under releases/), sorted.
.abc_remote_manifests <- function() {
    manifests <- character()
    token <- NULL

    repeat {
        page <- .abc_list_release_page(token)
        manifests <- c(manifests, paste0(page$prefixes, "manifest.json"))
        if (is.null(page$next_token)) {
            break
        }
        token <- page$next_token
    }

    sort(unique(manifests), method = "radix")
}

# One page of the S3 listing of release prefixes. `next_token` is NULL when
# the listing is complete.
.abc_list_release_page <- function(token = NULL) {
    url <- paste0(
        .abc_bucket_url, "/?list-type=2", "&prefix=releases/", "&delimiter=/"
    )
    if (!is.null(token)) {
        url <- paste0(url, "&continuation-token=", curl::curl_escape(token))
    }

    resp <- .abc_fetch_memory_retry(url)
    if (resp$status_code != 200L) {
        stop(
            "Could not list ABC Atlas releases (HTTP ", resp$status_code, ").",
            call. = FALSE
        )
    }

    doc <- xml2::xml_ns_strip(xml2::read_xml(resp$content))
    find_text <- function(xpath) {
        xml2::xml_text(xml2::xml_find_first(doc, xpath))
    }
    prefixes <- xml2::xml_text(
        xml2::xml_find_all(doc, "//CommonPrefixes/Prefix")
    )
    next_token <- if (identical(find_text("//IsTruncated"), "true")) {
        find_text("//NextContinuationToken")
    }

    list(prefixes = prefixes, next_token = next_token)
}


# Pick the manifest to load: the one last used with this cache if still valid,
# otherwise the latest available. Warns as abc_atlas_access does when that is
# not the newest release, or when it is newer than what was used before.
.abc_choose_manifest <- function(download_base, available) {
    latest <- utils::tail(available, 1L)
    last_used_file <- file.path(download_base, "_manifest_last_used.txt")

    if (file.exists(last_used_file)) {
        to_load <- paste(
            readLines(last_used_file, warn = FALSE), collapse = "\n"
        )
        if (!to_load %in% available) {
            warning("The manifest recorded as last used for this cache -- ",
                    to_load,
                " -- is not a valid manifest. Loading the latest version -- ",
                latest, " -- instead.",
                call. = FALSE
            )
            return(latest)
        }
        if (to_load != latest) {
            warning("You are loading ", to_load,
                ". A more up to date version of the ",
                "dataset -- ", latest, " -- exists. To use it, delete '",
                last_used_file, "' or point download_base to a new directory.",
                call. = FALSE
            )
        }
        return(to_load)
    }

    downloaded <- .abc_local_manifests(download_base)
    if (length(downloaded) > 0L && utils::tail(downloaded, 1L) != latest) {
        warning("You are loading ", latest,
            ", which is newer than the most recent ",
            "manifest previously used here (",
            utils::tail(downloaded, 1L), "). ",
            "Some data files may have changed between the two releases and ",
            "will be downloaded again (existing files are not overwritten).",
            call. = FALSE
        )
    }
    latest
}

.abc_check_manifest_version <- function(version) {
    rng <- .abc_manifest_compatibility
    if (is.null(version) || version < rng[1] || version >= rng[2]) {
        stop(
            "The manifest has version ", version, " but this package supports ",
            "versions on or after ", rng[1], " and before ", rng[2], "."
        )
    }
}

# ---- Files ------------------------------------------------------------------

# Resolve a file in the manifest. `file_name` is either a file name or
# "<name>/<kind>" for entries with several kinds (e.g. "WHB-10Xv3-Neurons/raw").
.abc_file_attributes <- function(cache, directory, file_name) {
    parts <- strsplit(file_name, "/", fixed = TRUE)[[1]]
    name <- parts[1]
    kind <- if (length(parts) > 1L) parts[2] else NULL

    dir_data <- cache$manifest$file_listing[[directory]]
    if (is.null(dir_data)) {
        stop(
            "Directory '", directory, "' not found in manifest ",
            cache$manifest_name, "."
        )
    }

    entry <- NULL
    for (sub_dir in names(dir_data)) {
        files_data <- dir_data[[sub_dir]][[name]]
        if (is.null(files_data)) next
        if (!is.null(files_data[["files"]])) {
            entry <- files_data[["files"]]
        } else if (!is.null(kind) && !is.null(files_data[[kind]])) {
            entry <- files_data[[kind]][["files"]]
        } else if (is.null(kind)) {
            stop(
                "File '", name, "' has several kinds in directory '", directory,
                "'. Specify the file name as one of: ",
                paste0(name, "/", names(files_data), collapse = ", ")
            )
        }
    }
    if (is.null(entry)) {
        stop("File '", name, "' not found in directory '", directory, "'.")
    }

    f <- entry[[1]]
    list(
        url = f$url, version = f$version, size = f$size,
        relative_path = f$relative_path, file_type = names(entry)[1],
        file_hash = f$file_hash,
        local_path = file.path(cache$download_base, f$relative_path)
    )
}

# TRUE when a metadata table is published for a manifest directory.
# This inspects the manifest only and never downloads the table.
.abc_metadata_exists <- function(cache, directory, file) {
    dir_data <- cache$manifest$file_listing[[directory]]
    if (is.null(dir_data)) {
        return(FALSE)
    }

    metadata <- dir_data[["metadata"]]
    !is.null(metadata) && !is.null(metadata[[file]])
}

# Expression matrix names and representations published for one manifest
# directory. This inspects the manifest only; it never downloads a file.
.abc_expression_inventory <- function(cache, directory) {
    dir_data <- cache$manifest$file_listing[[directory]]
    if (is.null(dir_data)) {
        stop(
            "Directory '", directory, "' not found in manifest ",
            cache$manifest_name, ".",
            call. = FALSE
        )
    }

    expression <- dir_data[["expression_matrices"]]
    if (is.null(expression) || !length(expression)) {
        return(data.frame(
            directory = character(),
            feature_matrix_label = character(),
            data_type = character(),
            stringsAsFactors = FALSE
        ))
    }

    out <- list()
    n <- 0L

    for (matrix_name in names(expression)) {
        entry <- expression[[matrix_name]]

        # Expression matrices normally have representation subentries such as
        # raw and log2. Only advertise representations that actually contain a
        # manifest file record.
        for (data_type in names(entry)) {
            if (identical(data_type, "files")) next

            representation <- entry[[data_type]]
            if (is.null(representation[["files"]])) next

            n <- n + 1L
            out[[n]] <- data.frame(
                directory = directory,
                feature_matrix_label = matrix_name,
                data_type = data_type,
                stringsAsFactors = FALSE
            )
        }
    }

    if (!length(out)) {
        return(data.frame(
            directory = character(),
            feature_matrix_label = character(),
            data_type = character(),
            stringsAsFactors = FALSE
        ))
    }

    do.call(rbind, out)
}

# Local path of a manifest file, downloading it first if needed (online only).
# A file on disk counts as present only if it is recorded as fully downloaded;
# otherwise it is fetched again, as abc_atlas_access does.
.abc_get_file_path <- function(
        cache, directory, file_name,
        force_download = FALSE,
        skip_hash_check = FALSE) {
    attrs <- .abc_file_attributes(cache, directory, file_name)
    path <- attrs$local_path
    if (dir.exists(path)) stop(path, " exists but is not a file.")

    if (cache$offline) {
        if (!file.exists(path)) {
            stop(
                "'", attrs$relative_path,
                "' is not available locally (offline = TRUE). ",
                "Run once with offline = FALSE to download it."
            )
        }
        .abc_report_cached(attrs)
        return(path)
    }

    base <- cache$download_base
    if (file.exists(path) && !.abc_is_recorded(base, path)) {
        force_download <- TRUE
    }
    if (force_download) {
        .abc_unrecord(base, path)
        unlink(path)
    }
    if (file.exists(path)) {
        .abc_report_cached(attrs)
        return(path)
    }

    message(
        "Downloading ", attrs$relative_path,
        " (", .abc_format_size(attrs$size), ")"
    )
    .abc_download(attrs$url, path, if (!skip_hash_check) attrs$file_hash)
    .abc_record(base, path)
    path
}

# Report that a requested file is served from the local cache. It mirrors the
# "Downloading ..." message so that reuse of earlier downloads is visible.
.abc_report_cached <- function(attrs) {
    message(
        "Using cached ", attrs$relative_path,
        " (", .abc_format_size(attrs$size), "); no download needed."
    )
}

# Download `url` to `dest` through a ".part" file that is resumed across retries
# and only moved into place once complete (and matching `md5`, if given).
.abc_download_progress <- function() {
    isTRUE(getOption("ABCatlas.download_progress", TRUE))
}
.abc_stop_if_download_interrupted <- function(
        result,
        part
) {
    if (
        length(result$success) != 1L ||
            !is.na(result$success)
    ) {
        return(invisible(FALSE))
    }

    partial_bytes <- if (file.exists(part)) {
        file.info(part)$size
    } else {
        0
    }

    partial_detail <- if (partial_bytes > 0) {
        sprintf(
            " (%s)",
            .abc_format_size(partial_bytes)
        )
    } else {
        ""
    }

    stop(
        "Download interrupted. Partial download retained for later resume",
        partial_detail,
        ".",
        call. = FALSE
    )
}

.abc_finalize_download_attempt <- function(
        result,
        part,
        dest,
        md5
) {
    ok <- isTRUE(result$success)

    if (
        ok &&
            !is.null(md5) &&
            unname(tools::md5sum(part)) != md5
    ) {
        message(
            "Checksum mismatch for ",
            basename(dest),
            "; downloading again."
        )

        unlink(part)
        return(FALSE)
    }

    if (!ok) {
        return(FALSE)
    }

    if (!file.rename(part, dest)) {
        stop(
            "Downloaded file could not be moved into place: ",
            dest,
            call. = FALSE
        )
    }

    TRUE
}

.abc_download_error_detail <- function(result) {
    if (
        !is.null(result$error) &&
            !is.na(result$error)
    ) {
        return(
            sprintf(
                ": %s",
                result$error
            )
        )
    }

    ""
}

.abc_report_download_retry <- function(
        attempt,
        result
) {
    message(
        "Download attempt ",
        attempt,
        " failed",
        .abc_download_error_detail(result),
        ". Retrying from the partial download."
    )

    invisible(NULL)
}

.abc_stop_download_failed <- function(
        url,
        max_tries,
        result
) {
    stop(
        "Could not download ",
        url,
        " in ",
        max_tries,
        " attempts",
        .abc_download_error_detail(result),
        ".",
        call. = FALSE
    )
}

.abc_download <- function(url, dest, md5 = NULL, max_tries = 10L) {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    part <- paste0(dest, ".part")

    for (i in seq_len(max_tries)) {
        # Abort a connection that cannot be opened within 10 s or that stays
        # below 1 kB/s for 60 s; the loop then retries and resumes the file.
        res <- curl::multi_download(
            url,
            part,
            resume = TRUE,
            progress = .abc_download_progress(),
            connecttimeout = 10,
            low_speed_limit = 1000,
            low_speed_time = 60
        )

        .abc_stop_if_download_interrupted(
            result = res,
            part = part
        )

        if (.abc_finalize_download_attempt(
            result = res,
            part = part,
            dest = dest,
            md5 = md5
        )) {
            return(invisible(dest))
        }

        if (i < max_tries) {
            .abc_report_download_retry(
                attempt = i,
                result = res
            )
        }
    }

    .abc_stop_download_failed(
        url = url,
        max_tries = max_tries,
        result = res
    )
}

.abc_format_size <- function(bytes) {
    format(structure(as.numeric(bytes), class = "object_size"),
        units = "auto",
        standard = "SI"
    )
}

.abc_available_space <- function(path = getwd()) {
    path <- path.expand(path)

    if (!dir.exists(path)) {
        stop(
            "Path does not exist: ",
            path,
            call. = FALSE
        )
    }

    path <- normalizePath(
        path,
        winslash = "/",
        mustWork = TRUE
    )

    if (.Platform$OS.type == "windows") {
        return(.abc_available_space_windows(path))
    }

    .abc_available_space_posix(path)
}

# Convert a string of decimal digits to bytes. Anything else (a header, a
# percentage, an empty string) gives NA instead of a coercion warning.
.abc_parse_bytes <- function(x, multiplier = 1) {
    x <- trimws(as.character(x))

    if (length(x) != 1L || is.na(x) || !grepl("^[0-9]+$", x)) {
        return(NA_real_)
    }

    as.numeric(x) * multiplier
}

# Run a system command and return its stdout lines, or character() if the
# command is missing, fails, or exits with a non-zero status.
.abc_command_output <- function(command, args) {
    tryCatch(
        system2(
            command,
            args,
            stdout = TRUE,
            stderr = FALSE
        ),
        warning = function(w) character(),
        error = function(e) character()
    )
}

# Free space, in bytes, on the filesystem holding `path` (macOS, Linux and
# other POSIX systems). `df -P` prints exactly one line per filesystem, so long
# device names (for example network mounts on HPC clusters) cannot wrap. The
# path is quoted because system2() passes arguments to the shell unquoted.
# Available space is the field immediately before the "NN%" capacity column,
# which tolerates spaces in the filesystem or mount-point names.
.abc_available_space_posix <- function(path) {
    out <- .abc_command_output(
        "df",
        c("-Pk", shQuote(path))
    )

    if (length(out) < 2L) {
        return(NA_real_)
    }

    fields <- strsplit(
        trimws(out[[length(out)]]),
        "[[:space:]]+"
    )[[1L]]

    capacity <- which(grepl("^[0-9]+%$", fields))

    if (!length(capacity) || capacity[[1L]] < 2L) {
        return(NA_real_)
    }

    .abc_parse_bytes(
        fields[[capacity[[1L]] - 1L]],
        multiplier = 1024
    )
}

.abc_available_space_windows <- function(path) {
    drive_letter <- substr(path, 1L, 1L)

    out <- .abc_command_output(
        "powershell",
        c(
            "-NoProfile",
            "-Command",
            sprintf("(Get-PSDrive -Name %s).Free", drive_letter)
        )
    )

    if (!length(out)) {
        return(NA_real_)
    }

    .abc_parse_bytes(out[[1L]])
}

# ---- Record of completed downloads (_downloaded_data.json) -------------------
# A JSON array of absolute paths, shared with abc_atlas_access.

.abc_record_path <- function(base) file.path(base, "_downloaded_data.json")

.abc_recorded <- function(base) {
    p <- .abc_record_path(base)
    if (!file.exists(p)) {
        return(character())
    }
    as.character(unlist(jsonlite::read_json(p)))
}

.abc_write_record <- function(base, paths) {
    jsonlite::write_json(unique(paths), .abc_record_path(base), pretty = TRUE)
}

.abc_is_recorded <- function(base, path) {
    normalizePath(path) %in% .abc_recorded(base)
}

.abc_record <- function(base, path) {
    .abc_write_record(base, c(.abc_recorded(base), normalizePath(path)))
}

.abc_unrecord <- function(base, path) {
    rec <- .abc_recorded(base)
    abs <- normalizePath(path, mustWork = FALSE)
    if (abs %in% rec) .abc_write_record(base, setdiff(rec, abs))
}

# ---- Metadata tables --------------------------------------------------------

# Read a metadata CSV into a data.frame. `col_classes` is a named vector of
# column classes (e.g. c(cell_label = "character")); other types are inferred.
# Integers too large for R's integer type are read as doubles; surrounding
# whitespace is kept, as pandas does.
.abc_read_csv <- function(path, col_classes = NULL, ...) {
    data.table::fread(
        path,
        colClasses = col_classes, na.strings = .abc_na_strings,
        integer64 = "double", strip.white = FALSE, data.table = FALSE,
        showProgress = FALSE, ...
    )
}
