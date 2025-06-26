setGeneric(
    "getGraph",
    signature = c("x"),
    function(x, ...) standardGeneric("getGraph")
)

setGeneric(
    "tableX",
    signature = c("x"),
    function(x, ...) standardGeneric("tableX")
)

setGeneric(
    "tableX<-",
    signature = c("x"),
    function(x, ..., value) standardGeneric("tableX<-")
)

setGeneric(
    "tableY",
    signature = c("x"),
    function(x, ...) standardGeneric("tableY")
)

setGeneric(
    "tableY<-",
    signature = c("x"),
    function(x, ..., value) standardGeneric("tableY<-")
)

setGeneric(
    "dictionary",
    signature = c("x"),
    function(x, ...) standardGeneric("dictionary")
)

setGeneric(
    "dictionary<-",
    signature = c("x"),
    function(x, ..., value) standardGeneric("dictionary<-")
)

setGeneric(
    "getEdgeList",
    signature = c("x"),
    function(x, ...) standardGeneric("getEdgeList")
)

setGeneric(
    "getFeaturePairs",
    signature = c("x"),
    function(x, ...) standardGeneric("getFeaturePairs")
)

#' @rdname weaveWeb
#' @export
#'
setGeneric(
    "weaveWeb",
    signature = c("x"),
    function(x, ...) standardGeneric("weaveWeb")
)

#' Run anansi
#'
#' @rdname getAnansi
#' @export
setGeneric(
    "anansi",
    signature = c("x"),
    function(x, ...) standardGeneric("anansi")
)

#' Bioc style plotting wrapper for anansi output
#'
#' @rdname plotAnansi
#' @export
setGeneric(
    "plotAnansi",
    signature = c("x"),
    function(x, ...) standardGeneric("plotAnansi")
)
