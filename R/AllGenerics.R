#' Get a listof edges
#' @name getEdgeList
#' @rdname getEdgeList
#' @param x input object
#' @param ... additional arguments
#' @export
getEdgeList        <- S7::new_generic("getEdgeList", "x")

#' Get a graph object.
#' @name getGraph
#' @seealso [getGraph.MultiFactor()]
#' @inheritParams getGraph.MultiFactor
#' @param x input
#' @export
getGraph           <- S7::new_generic("getGraph", "x")

#' Get a list of all pairs of features
#' @name getFeaturePairs
#' @rdname getFeaturePairs
#' @param x input object
#' @param ... additional arguments for specific methods
#' @export
getFeaturePairs    <- S7::new_generic("getFeaturePairs", "x")

#' Weave an AnansiWeb object
#' @name weaveWeb
#' @param x input object
#' @param ... additional arguments
#' @seealso [weaveWeb-methods()]
#' @export
weaveWeb           <- S7::new_generic("weaveWeb",   "x")

#' @name plotAnansi
#' @rdname plotAnansi
#' @export
#' @usage NULL
plotAnansi         <- S7::new_generic("plotAnansi", "x")

#' Apply a function on each pair of features
#' @name pairwiseApply
#' @rdname pairwiseApply
#' @param X input object
#' @param ... additional arguments
#' @export
pairwiseApply <- S7::new_generic("pairwiseApply", "X")

unfactor     <- S7::new_external_generic("S4Vectors", "unfactor", "x")
