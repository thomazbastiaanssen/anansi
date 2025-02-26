#' anansi wrapper for the MultiAssayExperiment class
#'
#' \code{getAnansi} provides a wrapper to execute the anansi pipeline on a
#' \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}
#' object. It applies the functions \code{\link{weaveWeb}} and
#' \code{\link{anansi}} in the given order.
#' @usage NULL
#'
#' @inheritParams getWeb
#'
#' @param ... additional parameters that can be passed to
#'   \code{\link{web}} or \code{\link{anansi}}.
#'
#' @details
#' This wrapper of \code{\link{anansi}} allows to perform a complete anansi
#' analysis directly on objects of class
#' \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}
#' . First, the assays specified by \code{assay.typeY} and \code{assay.typeX}
#' are passed to \code{\link{web}} to build an anansiWeb object.
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
#' # Load data
#' data("FMT_data", package = "anansi")
#'
#' # Convert to (Tree)SummarizedExperiment objects
#' metab_se <- SummarizedExperiment(assays = SimpleList(conc = as.matrix(FMT_metab)))
#' KO_tse <- TreeSummarizedExperiment(assays = SimpleList(counts = as.matrix(FMT_KOs)))
#'
#' # Select functions that are represented in the dictionary
#' keep <- row.names(KO_tse) %in% sort(unique(unlist(anansi_dic)))
#' KO_tse <- KO_tse[keep, ]
#'
#' # Remove features with less than 10% prevalence
#' KO_tse <- subsetByPrevalent(KO_tse,
#'   assay.type = "counts",
#'   prevalence = 0.1
#' )
#'
#' # Perform a centered log-ratio transformation on the functional counts assay
#' KO_tse <- transformAssay(KO_tse,
#'   assay.type = "counts",
#'   method = "clr",
#'   pseudocount = TRUE
#' )
#'
#' # Prepare colData
#' coldata <- FMT_metadata
#' rownames(coldata) <- coldata$Sample_ID
#' coldata <- coldata[match(colnames(KO_tse), rownames(coldata)), ]
#'
#' # Combine experiments into MultiAssayExperiment object
#' mae <- MultiAssayExperiment(
#'   experiments = ExperimentList(cpd = metab_se, ko = KO_tse),
#'   colData = coldata
#' )
#'
#' # Perform anansi analysis
#' out <- getAnansi(mae,
#'   experimentY = "cpd", experimentX = "ko",
#'   assay.typeY = "conc", assay.typeX = "clr",
#'   formula = ~Legend
#' )
#'
#' # Select significant interactions
#' out <- out[out$full_p.values < 0.05, ]
#'
#' # View subset of results
#' head(out, 5)
#'
#' @seealso
#' \code{\link{web}}
#' \code{\link{anansi}}
#'
#' @name getAnansi
#'
NULL

#' @rdname getAnansi
#' @export
setGeneric("getAnansi", signature = c("x"),
    function(x, ...) standardGeneric("getAnansi")
)

#' @rdname getAnansi
#' @export
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom SummarizedExperiment assay colData
setMethod("getAnansi",
  signature = c(x = "MultiAssayExperiment"),
  function(x, experimentY = 1, experimentX = 2, assay.typeY = "counts",
           assay.typeX = "counts", ...) {
    # Retrieve kwargs as list
    kwargs <- list(...)
    # Check fixed arguments
    fixed_args <- c("tableY", "tableX", "web", "metadata")
    remove <- names(kwargs) %in% fixed_args
    # If fixed arguments in kwargs, remove them
    if (any(remove)) {
      removed <- paste0(names(kwargs[remove]), sep = "'", collapse = ", '")
      kwargs <- kwargs[!remove]
      warning("The arguments '", removed, " should not be used, as they are ",
        "extracted from 'x'.",
        call. = FALSE
      )
    }
    # Check experiments
    mia:::.test_experiment_of_mae(x, experimentY)
    mia:::.test_experiment_of_mae(x, experimentX)
    # Extract experiments
    x1 <- x[[experimentY]]
    x2 <- x[[experimentX]]
    # Check assays
    mia:::.check_assay_present(assay.typeY, x1)
    mia:::.check_assay_present(assay.typeX, x2)
    # Extract assays
    t1 <- t(assay(x1, assay.typeY))
    t2 <- t(assay(x2, assay.typeX))
    # Extract colData (metadata)
    coldata <- colData(x)
    # Combine web.default args into list
    web_args <- c(list(x = experimentX, y = experimentY, tableX = t2, tableY = t1), kwargs)
    # Add kegg as a default link to match anansi()
    if(!"link" %in% names(web_args)) web_args[["link"]] <- kegg_link()
    keep <- names(web_args) %in% c("x", "y", "tableX", "tableY", "link")
    web_args <- web_args[keep]
    # Generate web object
    web <- do.call(web.default, web_args)
    # Combine anansi args into list
    anansi_args <- c(list(web = web, metadata = as.data.frame(coldata)), kwargs)
    keep <- names(anansi_args) %in% c(
      "web", "formula", "groups", "metadata", 
      "adjust.method", "verbose", "return.format"
    )
    anansi_args <- anansi_args[keep]
    # Generate anansi output
    out <- do.call(anansi, anansi_args)

    return(out)
  }
)
