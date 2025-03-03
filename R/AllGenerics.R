#' @rdname anansiLinkMap-methods
#' @usage NULL
#' @export
#' 
setGeneric("getGraph", signature = c("x"),
           function(x, ...) standardGeneric("getGraph")
)

#' weaveWeb wrapper for the MultiAssayExperiment class
#' @rdname getWeb
#' @usage NULL
#' @export
#' 
setGeneric("getWeb", signature = c("x"),
           function(x, ...) standardGeneric("getWeb")
)

#' anansi wrapper for the MultiAssayExperiment class
#'
#' @rdname getAnansi
#' @export
setGeneric("getAnansi", signature = c("x"),
           function(x, ...) standardGeneric("getAnansi")
)

#' miaViz style plotting wrapper for anansi output
#'
#' @rdname plotAnansi
#' @export
setGeneric("plotAnansi", signature = c("x"),
           function(x, ...) standardGeneric("plotAnansi")
)
