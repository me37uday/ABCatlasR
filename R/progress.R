#' @importFrom utils flush.console

.abc_progress_enabled <- function() {
  isTRUE(
    getOption(
      "ABCatlas.progress",
      interactive()
    )
  )
}

.abc_status_done <- function(text) {
  if (!.abc_progress_enabled()) {
    return(invisible(NULL))
  }

  cat(
    "\r",
    text,
    "    \n",
    sep = ""
  )

  flush.console()

  invisible(NULL)
}

.abc_confirm_download <- function() {
  repeat {
    answer <- trimws(
      readline("Proceed with download? [y/n]: ")
    )

    if (answer %in% c("y", "yes")) {
      return(TRUE)
    }

    if (answer %in% c("n", "no")) {
      cat(
        "Download cancelled. No expression files were downloaded.\n"
      )
      return(FALSE)
    }

    cat("Please enter y or n.\n")
  }
}

.abc_spinner_start <- function(label) {
  if (!.abc_progress_enabled()) {
    return(NULL)
  }

  rscript <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") {
      "Rscript.exe"
    } else {
      "Rscript"
    }
  )

  code <- paste(
    'label <- Sys.getenv("ABC_SPINNER_LABEL")',
    'i <- 1L',
    'repeat {',
    '  cat("\\r", label, strrep(".", i), "    ", sep = "")',
    '  flush.console()',
    '  Sys.sleep(1)',
    '  i <- if (i >= 4L) 1L else i + 1L',
    '}',
    sep = "; "
  )

  processx::process$new(
    command = rscript,
    args = c(
      "--vanilla",
      "-e",
      code
    ),
    stdout = "",
    stderr = NULL,
    env = c(
      ABC_SPINNER_LABEL = label
    ),
    cleanup = TRUE,
    windows_hide_window = TRUE
  )
}

.abc_spinner_stop <- function(spinner, text = NULL) {
  if (!is.null(spinner)) {
    try(
      spinner$kill(),
      silent = TRUE
    )
  }

  if (!.abc_progress_enabled()) {
    return(invisible(NULL))
  }

  cat(
    "\r",
    strrep(" ", 120L),
    "\r",
    sep = ""
  )

  if (!is.null(text)) {
    cat(
      text,
      "\n",
      sep = ""
    )
  }

  flush.console()

  invisible(NULL)
}
