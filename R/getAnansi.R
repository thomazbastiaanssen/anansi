#' anansi wrapper for the MultiAssayExperiment class
#'
#' \code{getAnansi} provides a wrapper to execute the anansi pipeline on a
#' \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}
#' object. It applies the functions \code{\link{weaveWebFromTables}} and
#' \code{\link{anansi}} in the given order.
#'
#' @param x a
#' \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}.
#'
#' @param experiment1 \code{Character scalar} or \code{numeric scalar}.
#'   Selects the experiment 1 from \code{experiments(x)} of
#'   \code{MultiassayExperiment} object. (Default: \code{1})
#'
#' @param experiment2 \code{Character scalar} or \code{numeric scalar}.
#'   Selects the experiment 2 from\code{experiments(x)} of
#'   \code{MultiAssayExperiment} object. (Default: \code{2})
#'
#' @param assay.type1 \code{Character scalar}. Specifies the name of the assay
#'   in experiment 1 to be used. (Default: \code{"counts"})
#'
#' @param assay.type2  \code{Character scalar}. Specifies the name of the
#'   assay in experiment 2 to be used. (Default: \code{"counts"})
#'
#' @param ... additional parameters that can be passed to
#'   \code{\link{weaveWebFromTables}} or \code{\link{anansi}}.
#'
#' @details
#' This wrapper of \code{\link{anansi}} allows to perform a complete anansi
#' analysis directly on objects of class
#' \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}
#' . First, the assays specified by \code{assay.type1} and \code{assay.type2}
#' are passed to \code{\link{weaveWebFromTables}} to build an anansiWeb object.
#' Next, this object is fed to the main \code{\link{anansi}} function to compute
#' interactions between the two assays.
#'
#' @return
#' If \code{return.format} is \code{"table"} (default), a wide format data.frame
#' intended to be compatible with `ggplot2`, or specialized plotting functions
#' (See \code{\link{[plotAnansi:plotAnansi]}}). If \code{return.format} is
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
#' data("dictionary", package = "anansi")
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
#'   experiments = ExperimentList(metabolites = metab_se, functions = KO_tse),
#'   colData = coldata
#' )
#'
#' # Perform anansi analysis
#' out <- getAnansi(mae,
#'   experiment1 = "metabolites", experiment2 = "functions",
#'   assay.type1 = "conc", assay.type2 = "clr",
#'   formula = ~Legend
#' )
#'
#' # Select significant interactions
#' out <- out[out$full_q.values < 0.1, ]
#'
#' # View subset of results
#' head(out, 5)
#'
#' @seealso
#' \code{\link{weaveWebFromTables}}
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
  function(x, experiment1 = 1, experiment2 = 2, assay.type1 = "counts",
           assay.type2 = "counts", ...) {
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
    mia:::.test_experiment_of_mae(x, experiment1)
    mia:::.test_experiment_of_mae(x, experiment2)
    # Extract experiments
    x1 <- x[[experiment1]]
    x2 <- x[[experiment2]]
    # Check assays
    mia:::.check_assay_present(assay.type1, x1)
    mia:::.check_assay_present(assay.type2, x2)
    # Extract assays
    t1 <- t(assay(x1, assay.type1))
    t2 <- t(assay(x2, assay.type2))
    # Extract colData (metadata)
    coldata <- colData(x)
    # Combine weaveWebFromTables args into list
    web_args <- c(list(tableY = t1, tableX = t2), kwargs)
    keep <- names(web_args) %in% c(
      "tableY", "tableX", "dictionary",
      "verbose", "mode", "prune", "max_sds"
    )
    web_args <- web_args[keep]
    # Generate web object
    web <- do.call(weaveWebFromTables, web_args)
    # Combine anansi args into list
    anansi_args <- c(list(web = web, metadata = as.data.frame(coldata)), kwargs)
    keep <- names(anansi_args) %in% c(
      "web", "formula", "groups", "metadata", "adjust.method",
      "verbose", "ignore_dictionary", "return.format"
    )
    anansi_args <- anansi_args[keep]
    # Generate anansi output
    out <- do.call(anansi, anansi_args)

    return(out)
  }
)
