#' anansi wrapper for the MultiAssayExperiment class
#'
#' \code{getAnansi} provides a wrapper to execute the anansi pipeline on a
#' \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}
#' object. It applies the functions \code{\link{weaveWeb}} and
#' \code{\link{anansi}} in the given order.
#' @inheritParams anansi
#' @inheritParams getWeb
#' @param ... additional parameters that can be passed to
#'   \code{\link{AnansiWeb}} or \code{\link{anansi}}.
#'
#' @details
#' This wrapper of \code{\link{anansi}} allows to perform a complete anansi
#' analysis directly on objects of class
#' \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}
#' . First, the assays specified by \code{assay.typeY} and \code{assay.typeX}
#' are passed to \code{\link{AnansiWeb}} to build an AnansiWeb object.
#' Next, this object is fed to the main \code{\link{anansi}} function to compute
#' interactions between the two assays.
#'
#' @return
#' If \code{return.format} is \code{"table"} (default), a wide format data.frame
#' intended to be compatible with `ggplot2`, or specialized plotting functions
#' (See \code{\link{plotAnansi}}). If \code{return.format} is
#' \code{"list"}, a list with aforementioned table, as well as input and
#' additional information. If \code{return.format} is \code{"raw"}, a list of
#' raw output (used for testing purposes).
#'
#' @examples
#'
#' # Import libraries
#' library(mia)
#' library(TreeSummarizedExperiment)
#' library(MultiAssayExperiment)
#'
#' web <- randomWeb(n_samples = 100)
#' mae <- as(web, "MultiAssayExperiment")
#'
#' # Perform anansi analysis
#' out <- getAnansi(mae,
#'   tableY = "y", tableX = "x",
#'   formula = ~ cat_ab
#' )
#'
#' # View subset of results
#' head(out, 5)
#'
#' @seealso
#' \code{\link{AnansiWeb}}
#' \code{\link{anansi}}
#'
#' @name getAnansi
#'
NULL

#' @rdname getAnansi
#' @export
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom SummarizedExperiment assay colData
#'
setMethod("getAnansi",
  signature = c(x = "MultiAssayExperiment"),
  function(x, tableY = 1, tableX = 2, formula, link = NULL, force_new = FALSE, ...) {
    # Retrieve kwargs as list
    kwargs <- list(...)
    # Check fixed arguments
    fixed_args <- c("web", "metadata")
    remove <- names(kwargs) %in% fixed_args
    # If fixed arguments in kwargs, remove them
    if (any(remove)) {
      removed <- paste0(names(kwargs[remove]), sep = "'", collapse = ", '")
      kwargs <- kwargs[!remove]
      stop("The arguments '", removed, " should not be used, as they are ",
        "extracted from 'x'.",
        call. = FALSE
      )
    }
    # Generate web object
    w <- getWeb(x, tableY, tableX, link, ...)

    # Generate anansi output
    out <- anansi(web = w, formula = formula, ...)

    return(out)
  }
)
