#' @rdname getWeb
#' @export
#'
#' @param x a \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}.
#'
#' @param tableY,tableX
#' \code{Character scalar} or \code{numeric scalar}. Selects experiment
#' corresponding to \code{tableY} and \code{tableX} from \code{experiments(x)}
#' of \code{MultiAssayExperiment} object by name or index, name is recommended.
#' (Default slots: \code{Y = 1}, \code{X = 2}).
#'
#' @param link One of the following:
#' \itemize{
#'  \item Character scalar with value "none"
#'  \item data.frame with two columns
#'  \item list with two such data.frames
#' }
#'
#' @param force_new \code{boolean} If x already has a dictionary \code{Matrix}
#' in metadata, ignore it and generate a new object anyway? (Default: FALSE).
#' @param ... additional parameters passed to \code{\link{AnansiWeb}}.
#' @param tableY,tableX
#' \code{Character scalar} or \code{numeric scalar}. Selects experiment
#' corresponding to \code{tableY} and \code{tableX} from \code{experiments(x)}
#' of \code{MultiAssayExperiment} object by name or index, name is recommended.
#' (Default slots: \code{Y = 1}, \code{X = 2}).
#'
#' @returns an \code{AnansiWeb} object, with sparse binary biadjacency matrix
#' with features from \code{y} as rows and features from \code{x} as columns in
#' \code{dictionary} slot. If x already contains a dictionary in metadata, use
#' that one, unless \code{force_new = TRUE}.
#'
#' @importFrom MultiAssayExperiment MultiAssayExperiment metadata metadata<-
#' @importFrom SummarizedExperiment assay colData
#' @importClassesFrom Matrix Matrix
#'
#' @details
#' This wrapper of \code{\link{weaveWeb}} allows to generate an
#' \code{AnansiWeb} S4 object directly from objects of class
#' \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}
#' . First, the assays specified by \code{assay.typeY} and \code{assay.typeX}
#' are passed to \code{\link{AnansiWeb}} to build an AnansiWeb object.
#'
#' @examples
#' # Make a random anansiWeb
#' web <- randomWeb()
#'
#' # Combine experiments into MultiAssayExperiment object
#' mae <- as(web, "MultiAssayExperiment")
#'
#' # Back to AnansiWeb
#' outWeb <- getWeb(mae, tableY = "y", tableX = "x")
#'
setMethod("getWeb", signature = c(x = "MultiAssayExperiment"), function(
        x, tableY = 1, tableX = 2, link = NULL, force_new = FALSE, ...
        ) {

    # Check experiments
    mia:::.test_experiment_of_mae(x, tableY)
    mia:::.test_experiment_of_mae(x, tableX)
    y_id <- names(experiments(mae)[tableY])
    x_id <- names(experiments(mae)[tableX])

    # Extract assays
    tableY <- t(assay(x, y_id))
    tableX <- t(assay(x, x_id))

    if(!force_new){
    # Check if x already contains a dictionary
    m <- metadata(x)

    if(is.null(link))
        d <- "dictionary"
    if(length(link) == 1L && (is.character(link) || is.numeric(link)) )
        d <- link
    if(d %in% names(m))
        return(AnansiWeb( tableX, tableY, m[[d]]))
    }
    # Generate web object
    weaveWeb.default(x, y, link, tableX, tableY, metadata = colData(x), ...)
}
)
