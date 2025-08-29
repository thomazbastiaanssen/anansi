#' Get a graph object out of an MultiFactor.
#' @name getGraph
#' @rdname getGraph
#' @aliases getGraph
#' @param x `MultiFactor`
#' @param format
#' `Character scalar`, controls output format by package name.
#' `"igraph"` and `"graph"` are supported.
#' @param ... additional arguments (currently not used).
#' @importFrom igraph graph_from_data_frame as_graphnel
#' @returns a specified graph object.
#' @seealso [igraph::graph_from_data_frame()] and
#'     [igraph::as_graphnel()], which are used under the hood, from
#'     [igraph::igraph()] package.
#' @export
#' @examples
#' # Generate a regular igraph object
#' g <- getGraph(kegg_link())
#' plot(g)
#'
#'
#' # Output formats
#' ec2cpd <- kegg_link()[["ec2cpd"]]
#'
#' getGraph(ec2cpd, format = "graph")
#' getGraph(ec2cpd, format = "igraph")
#'
S7::method(getGraph, MultiFactor) <- function(x, format = "igraph") {

        g <- graph_from_data_frame(getEdgeList(x), directed = FALSE)

        switch(
            format,
            "igraph" = {
            },
            "graph" = g <- as_graphnel(g)
        )
        return(g)
    }

#' @export
#' @name getGraph
#' @rdname getGraph
#'
S7::method(getGraph, S7::class_list | S7::class_data.frame) <- function(x, format = "igraph") {
    getGraph(asMultiFactor(x), format)
}

