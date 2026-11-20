#' Find the folder where ABC Atlas files are stored
#'
#' Returns the folder where ABC Atlas files are downloaded and kept for
#' reuse. It is the first of these that is set:
#' \enumerate{
#'   \item the explicit \code{path} argument;
#'   \item \code{getOption("ABCatlasR.download_base")};
#'   \item the \code{ABC_ATLAS_CACHE} environment variable;
#'   \item a per-user cache directory,
#'     \code{tools::R_user_dir("ABCatlasR", "cache")}.
#' }
#' The folder is created if it does not exist yet.
#'
#' @param path Folder to use. If \code{NULL} or empty, the other settings
#'   above are checked in order.
#' @return The full path of the folder.
#' @examples
#' # Use a temporary cache directory.
#' cache_dir <- file.path(tempdir(), "ABCatlasR-example")
#' abc_download_dir(cache_dir)
#' @export
abc_download_dir <- function(path = NULL) {
    is_set <- function(x) {
        !is.null(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    }

    dir <- NULL
    if (is_set(path)) {
        dir <- path
    } else if (is_set(getOption("ABCatlasR.download_base"))) {
        dir <- getOption("ABCatlasR.download_base")
    } else {
        env <- Sys.getenv("ABC_ATLAS_CACHE", unset = "")
        if (is_set(env)) {
            dir <- env
        } else {
            dir <- tools::R_user_dir("ABCatlasR", which = "cache")
        }
    }

    dir <- path.expand(dir)
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    normalizePath(dir, mustWork = FALSE)
}
