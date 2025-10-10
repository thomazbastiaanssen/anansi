#' @title Extract a graph object from a MultiFactor.
#' @name getGraph.anansi::MultiFactor
#' @rdname MultiFactor-methods
#' @param x `MultiFactor`
#' @param format
#' `Character scalar`, controls output format by package name.
#' `"igraph"` and `"graph"` are supported.
#' @importFrom igraph graph_from_data_frame as_graphnel
#' @returns a specified graph object.
#' @seealso [igraph::graph_from_data_frame()] and
#'     [igraph::as_graphnel()], which are used under the hood, from
#'     [igraph::igraph()] package.
#' @examples
#' # Generate an igraph object
#' g <- getGraph(x = kegg_link(), format = "igraph")
#' plot(g)
#'
NULL

#' @export
S7::method(getGraph, MultiFactor) <- function(
        x, format = "igraph"
) `getGraph.anansi::MultiFactor`(x, format)

`getGraph.anansi::MultiFactor` <- function(x, format = "igraph") {
    g <- graph_from_data_frame(getEdgeList(x), directed = FALSE)

    switch(format,
           "igraph" = {
           },
           "graph" = g <- as_graphnel(g)
    )
    return(g)
}

#' @method getGraph list
#' @name getGraph.list
#' @rdname getGraph
#'
S7::method(getGraph, S7::class_list | S7::class_data.frame) <- function(
        x, format = "igraph"
) {
    getGraph(asMultiFactor(x), format)
}
