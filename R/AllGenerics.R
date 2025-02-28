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