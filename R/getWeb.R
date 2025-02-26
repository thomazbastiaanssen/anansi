#' weaveWeb wrapper for the MultiAssayExperiment class
#' @rdname getWeb
#' @usage NULL
#' @export
#' 
setGeneric("getWeb", signature = c("x"),
           function(x, ...) standardGeneric("getWeb")
)

#' @rdname getWeb
#' @export
#' 
#' @param x a \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}.
#' 
#' @param experimentY,experimentX 
#' \code{Character scalar} or \code{numeric scalar}. Selects experiment 
#' corresponding to \code{tableY} and \code{tableX} from \code{experiments(x)} 
#' of \code{MultiAssayExperiment} object by name or index, name is recommended. 
#' (Default slots: \code{Y = 1}, \code{X = 2}).
#' 
#' @param assay.typeY,assay.typeX 
#' \code{Character scalar}. Specifies the name of the assay in experiments Y and
#' X to be used. (Default: \code{"counts"})
#' 
#' @param link One of the following:
#' \itemize{
#'  \item Character scalar with value "none"
#'  \item data.frame with two columns
#'  \item list with two such data.frames
#' }
#' 
#' @param ... additional parameters that can be passed to \code{\link{web}}.
#'
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom SummarizedExperiment assay colData
#' 
#' @details
#' This wrapper of \code{\link{weaveWeb}} allows to generate an 
#' \code{anansiWeb} S4 object directly from objects of class
#' \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}
#' . First, the assays specified by \code{assay.typeY} and \code{assay.typeX}
#' are passed to \code{\link{web}} to build an anansiWeb object.
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
#' outWeb <- getWeb(mae,
#'   experimentY = "cpd", experimentX = "ko",
#'   assay.typeY = "conc", assay.typeX = "clr"
#' )
#' 
setMethod("getWeb",
          signature = c(x = "MultiAssayExperiment"),
          function(x, experimentY = 1, experimentX = 2, 
                   assay.typeY = "counts", assay.typeX = "counts", 
                   link = NULL, ...) {
            # Retrieve kwargs as list
            kwargs <- list(...)
            # Check experiments
            mia:::.test_experiment_of_mae(x, experimentY)
            mia:::.test_experiment_of_mae(x, experimentX)
            # Extract experiments
            tableY <- x[[experimentY]]
            tableX <- x[[experimentX]]
            # Check assays
            mia:::.check_assay_present(assay.typeY, tableY)
            mia:::.check_assay_present(assay.typeX, tableX)
            # Extract assays
            tableY <- t(assay(tableY, assay.typeY))
            tableX <- t(assay(tableX, assay.typeX))
            # Combine web.default args into list
            web_args <- c(list(x = experimentX, y = experimentY, 
                               tableX = tableX, tableY = tableY), kwargs)
            # Add kegg as a default link to match anansi()
            if(!"link" %in% names(web_args)) web_args[["link"]] <- kegg_link()
            keep <- names(web_args) %in% c("x", "y", "tableX", "tableY", "link")
            web_args <- web_args[keep]
            # Generate web object
            do.call(web.default, web_args)
          }
)
