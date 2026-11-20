#' ABCatlasR: access the Allen Brain Cell Atlas from R
#'
#' ABCatlasR reads cell metadata and expression of the Allen Brain Cell (ABC)
#' Atlas directly in R, from a fixed ABC Atlas release, and returns them as a
#' sparse matrix bundle, a SingleCellExperiment or a SpatialExperiment.
#'
#' The workflow has four steps:
#' \enumerate{
#'   \item [list_datasets()] lists the supported dataset routes.
#'   \item [load_data()] loads the cell and gene metadata of one route.
#'   \item [plan_fetch()] optionally reports the expression files a
#'     selection of cells needs, their sizes and the disk space required.
#'   \item [fetch_data()] downloads and reads the expression of the
#'     selected cells and genes.
#' }
#'
#' Downloaded files are cached and reused; [abc_download_dir()] shows and
#' sets the cache location, and [abc_cache_connect()] gives direct access
#' to the files of a release. The vignette walks through a complete example:
#' \code{vignette("ABCatlasR", package = "ABCatlasR")}.
#'
#' @keywords internal
"_PACKAGE"
