#' Get a graph object out of an AnansiLinkMap.
#' @rdname getGraph
#' @param x `AnansiLinkMap`
#' @param format
#' `Character scalar`, controls output format by package name.
#' `"igraph"` and `"graph"` are supported.
#'
#' @param ... additional arguments (currently not used).
#' @importFrom igraph graph_from_data_frame as_graphnel
#' @seealso
#' [igraph::graph_from_data_frame()] and
#' [igraph::as_graphnel()], which are used under the hood,
#' from [igraph::igraph()] package.
#' @export
#' @examples
#' # Generate a regular igraph object
#' g <- getGraph( kegg_link() )
#' plot(g)
#'
#' # Output formats
#' getGraph( ec2cpd, format =  "graph" )
#' getGraph( ec2cpd, format = "igraph" )
#'
#'
setMethod("getGraph", "AnansiLinkMap",
          function(x, format = "igraph", ...) {
            validObject(x)

            g <- graph_from_data_frame(getEdgeList(x), directed = FALSE)

            switch(format,
                   "igraph" = {},
                   "graph"  = g <- as_graphnel(g)
            )
            return(g)
          }
)

#' @export
#' @rdname getGraph
#'
setMethod("getGraph", "list",
          function(x, format = "igraph", ...) getGraph( asLinkMap(x), format)
)
