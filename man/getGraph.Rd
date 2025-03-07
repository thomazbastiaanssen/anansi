% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods-getGraph.R
\name{getGraph}
\alias{getGraph}
\alias{getGraph,AnansiLinkMap-method}
\alias{getGraph,list-method}
\title{extract or make a graph}
\usage{
\S4method{getGraph}{AnansiLinkMap}(x, format = "igraph", ...)

\S4method{getGraph}{list}(x, format = "igraph", ...)
}
\arguments{
\item{x}{\code{AnansiLinkMap}}

\item{...}{additional arguments (currently not used).}

\item{format}{\verb{Character scalar}, controls output format by package name.
\code{"igraph"} and \code{"graph"} are supported.}
}
\description{
extract or make a graph

Get a graph object out of an AnansiLinkMap.
}
\examples{
# Generate a regular igraph object
g <- getGraph( kegg_link() )
plot(g)

# Output formats
getGraph( ec2cpd, format =  "graph" )
getGraph( ec2cpd, format = "igraph" )


}
\seealso{
\code{\link[igraph:graph_from_data_frame]{igraph::graph_from_data_frame()}} and
\code{\link[igraph:as_graphnel]{igraph::as_graphnel()}}, which are used under the hood,
from \code{\link[igraph:aaa-igraph-package]{igraph::igraph()}} package.
}
