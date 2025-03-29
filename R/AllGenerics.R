setGeneric("getGraph",
    signature = c("x"),
    function(x, ...) standardGeneric("getGraph")
)

setGeneric("tableX",
    signature = c("x"),
    function(x, ...) standardGeneric("tableX")
)

setGeneric("tableX<-",
    signature = c("x"),
    function(x, ..., value) standardGeneric("tableX<-")
)

setGeneric("tableY",
    signature = c("x"),
    function(x, ...) standardGeneric("tableY")
)

setGeneric("tableY<-",
    signature = c("x"),
    function(x, ..., value) standardGeneric("tableY<-")
)

setGeneric("dictionary",
    signature = c("x"),
    function(x, ...) standardGeneric("dictionary")
)

setGeneric("dictionary<-",
    signature = c("x"),
    function(x, ..., value) standardGeneric("dictionary<-")
)

setGeneric("getEdgeList",
    signature = c("x"),
    function(x, ...) standardGeneric("getEdgeList")
)

#' weaveWeb wrapper for the MultiAssayExperiment class
#' @rdname getWeb
#' @usage NULL
#' @export
#'
setGeneric("getWeb",
    signature = c("x"),
    function(x, ...) standardGeneric("getWeb")
)

#' anansi wrapper for the MultiAssayExperiment class
#'
#' @rdname getAnansi
#' @export
setGeneric("getAnansi",
    signature = c("x"),
    function(x, ...) standardGeneric("getAnansi")
)

#' miaViz style plotting wrapper for anansi output
#'
#' @rdname plotAnansi
#' @export
setGeneric("plotAnansi",
    signature = c("x"),
    function(x, ...) standardGeneric("plotAnansi")
)
