#' Get a listof edges
#' @name getEdgeList
#' @rdname getEdgeList
#' @param x input object
#' @param ... additional arguments
#' @examples
#' x <- randomMultiFactor(n_features = 10)
#' getEdgeList(x)
#' @returns a two-column data.frame that lists the content of each entry in the
#' input MultiFactor
#' @export
getEdgeList        <- S7::new_generic("getEdgeList", "x")

#' Get a graph object.
#' @name getGraph
#' @seealso [getGraph.MultiFactor()]
#' @inheritParams getGraph.MultiFactor
#' @param x input
#' @examples
#' x <- randomMultiFactor(n_features = 10)
#' getGraph(x)
#' @returns a specified graph object.
#' @export
getGraph           <- S7::new_generic("getGraph", "x")

#' Get a list of all pairs of features
#' @name getFeaturePairs
#' @rdname getFeaturePairs
#' @param x input object
#' @param ... additional arguments for specific methods
#' @returns an list of two-column data.frames that represent all feature pairs.
#' @examples
#' x <- randomWeb(10)
#' getFeaturePairs(x)
#'
#' @export
getFeaturePairs    <- S7::new_generic("getFeaturePairs", "x")

#' Weave an AnansiWeb object
#' @name weaveWeb
#' @param x input object
#' @param ... additional arguments
#' @seealso [weaveWeb-methods()]
#' @returns an AnansiWeb object
#' @examples
#' # Setup demo tables
#' ec2ko <- kegg_link()[["ec2ko"]]
#' ec2cpd <- kegg_link()[["ec2cpd"]]
#'
#' # Basic usage
#' weaveWeb(cpd ~ ko, link = kegg_link())
#' weaveWeb(x = "ko", y = "ec", link = ec2ko)
#' weaveWeb(ec ~ cpd, link = ec2cpd)
#'
#' # A wrapper is available for kegg ko, ec and cpd data
#' generic <- weaveWeb(cpd ~ ko, link = kegg_link())
#' kegg_wrapper <- weaveKEGG(cpd ~ ko)
#'
#' identical(generic, kegg_wrapper)
#'
#' # The following are equivalent to transposition:
#' a <- weaveWeb(ko ~ cpd, link = kegg_link())@dictionary
#' b <- weaveWeb(cpd ~ ko, link = kegg_link())@dictionary
#'
#' identical(a, Matrix::t(b))
#'
#' @export
weaveWeb           <- S7::new_generic("weaveWeb",   "x")

#' @name plotAnansi
#' @rdname plotAnansi
#' @export
#' @returns a figure that can be further modified using the `ggplot2` suite
#' @usage NULL
plotAnansi         <- S7::new_generic("plotAnansi", "x")

#' Apply a function on each pair of features
#' @name pairwiseApply
#' @rdname pairwiseApply
#' @param X input object
#' @param ... additional arguments
#' @returns
#' a list containing the output of applying the function to each feature pair.
#' See `?base::mapply()`
#' @examples
#' web <- randomWeb(10)
#'
#' pairwiseApply(
#'     X = web,
#'     FUN = function(x, y) cor(x, y),
#'     MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE
#' )
#' @export
pairwiseApply <- S7::new_generic("pairwiseApply", "X")

