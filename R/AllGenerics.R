#' Get a listof edges
#' @name getEdgeList
#' @rdname getEdgeList
#' @param x input object
#' @param ... additional arguments
#' @examples
#' x <- randomMultiFactor(n_features = 10)
#' getEdgeList(x)
#' @returns a two-column data.frame that lists the content of each entry in the
#'     input MultiFactor
#' @export
getEdgeList <- S7::new_generic("getEdgeList", "x")

#' Get a graph object.
#' @name getGraph
#' @seealso [`MultiFactor-methods`]
#' @param x input
#' @param ... additional arguments (currently not used).
#' @examples
#' # Show methods
#' getGraph
#' @returns a specified graph object.
#' @export
getGraph <- S7::new_generic("getGraph", "x")

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
getFeaturePairs <- S7::new_generic("getFeaturePairs", "x")

#' Weave an AnansiWeb object
#' @name weaveWeb
#' @rdname weaveWeb-generic
#' @param x input object
#' @param ... additional arguments
#' @seealso [weaveWeb-methods()]
#' @returns an `AnansiWeb` object, with sparse binary biadjacency matrix
#'     with features from `y` as rows and features from `x` as columns in
#'     `dictionary` slot.
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
#' a <- weaveWeb(ko ~ cpd, link = kegg_link()) |> dictionary()
#' b <- weaveWeb(cpd ~ ko, link = kegg_link()) |> dictionary()
#'
#' identical(a, Matrix::t(b))
#'
#' @export
weaveWeb <- S7::new_generic("weaveWeb", "x")

#' @name plotAnansi
#' @rdname plotAnansi
#' @aliases plotAnansi-generic
#' @export
#' @returns a figure that can be further modified using the `ggplot2` suite
#' @usage NULL
plotAnansi <- S7::new_generic("plotAnansi", "x")

#' Apply a function on each pair of features
#' @name pairwiseApply
#' @rdname pairwiseApply
#' @aliases pairwiseApply-generic
#' @param X input object
#' @param ... additional arguments
#' @returns
#' a list containing the output of applying the function to each feature pair.
#'     See `?base::mapply()`
#' @examples
#' web <- randomWeb(10)
#'
#' # For each feature pair, was the value for x higher than the value for y?
#' pairwiseApply(
#'     X = web,
#'     FUN = function(x, y) x > y,
#'     MoreArgs = NULL, SIMPLIFY = FALSE, USE.NAMES = TRUE
#'     )
#'
#' # Run cor.test() on each pair of features
#' pairwiseApply(
#'     X = web,
#'     FUN = function(x, y) cor.test(x, y),
#'     MoreArgs = NULL, SIMPLIFY = FALSE, USE.NAMES = TRUE
#' )
#'
#' @export
pairwiseApply <- S7::new_generic("pairwiseApply", "X")

#' @export
unfactor <- S7::new_external_generic(
    package = "S4Vectors", name = "unfactor", "x"
    )

#' @export
metadata <- S7::new_external_generic(
    package = "S4Vectors", name = "metadata", "x"
    )

#' @export
`metadata<-` <- S7::new_external_generic(
    package = "S4Vectors", name = "metadata<-", "x"
    )


#' Get tableX
#' @name tableX
#' @rdname tableX
#' @aliases tableX-generic
#' @param x input object
#' @param ... additional arguments
#' @returns `tableX` slot of x
#' @examples
#' x <- randomWeb(10)
#' tableX(x)
#' @export
#'
tableX <- S7::new_generic("tableX", "x")

#' Get tableY
#' @name tableY
#' @aliases tableY-generic
#' @rdname tableY
#' @param x input object
#' @param ... additional arguments
#' @returns `tableY` slot of x
#' @examples
#' x <- randomWeb(10)
#' tableY(x)
#' @export
#'
tableY <- S7::new_generic("tableY", "x")

#' Get dictionary
#' @name dictionary
#' @rdname dictionary
#' @aliases dictionary-generic
#' @param x input object
#' @param ... additional arguments
#' @returns `dictionary` slot of x
#' @examples
#' x <- randomWeb(10)
#' dictionary(x)
#'
#' @export
#'
dictionary <- S7::new_generic("dictionary", "x")
