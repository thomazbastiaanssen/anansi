#' @name generic.getEdgeList
#' @rdname MultiFactor
#' @usage NULL
#' @export
getEdgeList        <- S7::new_generic("getEdgeList", "x")

#' @name generic.getGraph
#' @rdname getGraph
#' @usage NULL
#' @export
getGraph           <- S7::new_generic("getGraph", "x")

#' @export
#' @rdname AnansiWeb
#' @usage NULL
#' @name generic.getFeaturePairs
getFeaturePairs    <- S7::new_generic("getFeaturePairs", "x")

#' @export
#' @rdname weaveWeb
#' @usage NULL
#' @name generic.weaveWeb
weaveWeb           <- S7::new_generic("weaveWeb",   "x")

#' @export
#' @rdname plotAnansi
#' @name generic.plotAnansi
#' @usage NULL
plotAnansi         <- S7::new_generic("plotAnansi", "x")

#' @export
#' @rdname AnansiWeb
#' @name generic.pairwiseApply
#' @usage NULL
pairwiseApply <- S7::new_generic("pairwiseApply", "X")

unfactor     <- S7::new_external_generic("S4Vectors", "unfactor", "x")
