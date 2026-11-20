#' ABCatlasR: access the Allen Brain Cell Atlas from R
#'
#' ABCatlasR reads cell metadata and gene expression from the Allen Brain
#' Cell (ABC) Atlas into R. Data come from a fixed ABC Atlas release, so
#' results can be reproduced, and expression is returned as a sparse matrix
#' with its cell and gene tables, a SingleCellExperiment or a
#' SpatialExperiment.
#'
#' The workflow has four steps:
#' \enumerate{
#'   \item [list_datasets()] lists the datasets you can load.
#'   \item [load_data()] loads the cell and gene metadata of one dataset.
#'   \item [plan_fetch()] optionally reports the expression files a
#'     selection of cells needs, their sizes and the disk space required.
#'   \item [fetch_data()] downloads and reads the expression of the
#'     selected cells and genes.
#' }
#'
#' Downloaded files are kept in a local folder and reused.
#' [abc_download_dir()] shows and sets that folder, and
#' [abc_cache_connect()] gives direct access to the files of a release.
#' The vignette walks through a complete example:
#' \code{vignette("ABCatlasR", package = "ABCatlasR")}.
#'
#' @author
#' In the order of the package \code{DESCRIPTION}, which also lists their
#' e-mail addresses and ORCID identifiers:
#' \itemize{
#'   \item Uday Rangaswamy
#'   \item Mahdi Komeili (maintainer, \email{mahdi.komeili1999@gmail.com})
#'   \item Federico Caretti
#'   \item Idris Kouadri Boudjelthia
#'   \item Maria Musarra Pizzo
#'   \item Amir Fardin
#'   \item Andres Hernandez-Clavijo
#'   \item Remo Sanges
#' }
#'
#' @keywords internal
"_PACKAGE"
