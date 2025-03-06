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
#' @importFrom MultiAssayExperiment MultiAssayExperiment metadata experiments
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
    y_id <- names(experiments(x)[tableY])
    x_id <- names(experiments(x)[tableX])

    # Extract assays
    tY <- t(assay(x, y_id))
    tX <- t(assay(x, x_id))

    if(!force_new){
    # Check if x already contains a dictionary
    m <- metadata(x)

    if(is.null(link))
        d <- "dictionary"
    if(valid_selection(link, m))
        d <- link
    if(d %in% names(m))
        return(AnansiWeb(tableX = tX, tableY = tY, dictionary = m[[d]],
                         metadata = colData(x), ...) )
    }
    # Generate web object
    weaveWeb.default(x = x_id, y = y_id, link = link,
                     tableX = tX, tableY = tY, metadata = colData(x), ...)
    }
)

#' TRUE if i can select in x
#' @noRd
#' @param i \code{Character or numeric scalar}. Index to check.
#' @param x object to check i in
#' @returns TRUE if i selects in x, FALSE otherwise
#'
valid_selection <- function(i, x) {
    # Need to be length 1.
    if(length(i) != 1L) FALSE

    # If numeric, needs to be within element length of x
    if(is.numeric(i)) i <= length(x)
    # If character, x needs to be named and i needs to be within those names
    if(is.character(i)){
        if(is.null(names(x))) FALSE

        !is.na(match(i, names(x)))
    }
    # If that didn't work, invalid selection. return FALSE.
    FALSE
}
