#' Get a graph object out of an anansiLinkMap. 
#' @rdname getGraph
#' @param x \code{anansiLinkMap}
#' @param format 
#' \code{Character scalar}, controls output format by package name. 
#' \code{"igraph"} and \code{"graph"} are supported. 
#' 
#' @param ... additional arguments (currently not used). 
#' @importFrom igraph graph_from_data_frame as_graphnel
#' @seealso 
#' \code{\link[igraph:graph_from_data_frame]{graph_from_data_frame}} and 
#' \code{\link[igraph:as_graphnel]{as_graphnel}}, which are used under the hood, 
#' from \code{\link[igraph:igraph]{igraph}} package.
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
setMethod("getGraph", "anansiLinkMap", 
          function(x, format = "igraph", ...) {
            validObject(x)
            
            g <- graph_from_data_frame(x@edgelist, directed = FALSE)
            
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
