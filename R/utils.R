#' Resolve the ABC Atlas download/cache directory
#'
#' Determines where ABC Atlas data is downloaded and cached. The location is
#' resolved with the following precedence (first non-empty wins):
#' \enumerate{
#'   \item the explicit \code{path} argument;
#'   \item \code{getOption("ABCatlas.download_base")};
#'   \item the \code{ABC_ATLAS_CACHE} environment variable;
#'   \item a per-user cache directory,
#'     \code{tools::R_user_dir("ABCatlas", "cache")}.
#' }
#' The resolved directory is created if it does not yet exist and returned as an
#' absolute path. This lets users point downloads at a custom folder without
#' depending on the working directory.
#'
#' @param path Optional explicit directory. If \code{NULL} or empty, the other
#'   sources above are consulted in order.
#' @return An absolute path to the (existing) download directory.
#' @examples
#' # Use a temporary cache directory.
#' cache_dir <- file.path(tempdir(), "ABCatlas-example")
#' abc_download_dir(cache_dir)
#' @export
abc_download_dir <- function(path = NULL) {
  is_set <- function(x) !is.null(x) && length(x) == 1L && !is.na(x) && nzchar(x)

  dir <- NULL
  if (is_set(path)) {
    dir <- path
  } else if (is_set(getOption("ABCatlas.download_base"))) {
    dir <- getOption("ABCatlas.download_base")
  } else {
    env <- Sys.getenv("ABC_ATLAS_CACHE", unset = "")
    if (is_set(env)) {
      dir <- env
    } else {
      dir <- tools::R_user_dir("ABCatlas", which = "cache")
    }
  }

  dir <- path.expand(dir)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(dir, mustWork = FALSE)
}
