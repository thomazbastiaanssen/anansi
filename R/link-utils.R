#' make an igraph mapping across the content of link
#' @details Make an \code{igraph} graph from the  content of \code{link}.  
#' @param link \code{data.frame} with exactly two, named, columns, or a 
#' @returns an \code{igraph} graph. 
#' @importFrom igraph graph_from_edgelist
#' @seealso \code{\link[igraph:igraph]{igraph}}
#' @export
#' @examples 
#' graph_from_link( ec2cpd )
#' 
#' g <- graph_from_link( kegg_link() )
#' plot(g)
#' 
graph_from_link <- function(link) {
  stopifnot("'link' must be one or several data.frame objects." =
              validLink(link))
  linknames <- colnames_link(link)
  m <- t(as.data.frame(linknames))
  
  graph_from_edgelist(m, directed = FALSE)
  
}

#' Is this one or several data.frames with exactly two columns that are named?
#' @noRd
validLink <- function(link) validLinkDF(link) || 
  all(unlist(lapply(link, validLinkDF)))

#' Is this a data.frame with exactly two columns that are named?
#' @noRd
validLinkDF <- function(x) is.data.frame(x) && 
  NCOL(x) == 2L &&
  length(colnames(x)) == 2L

#' get colnames
#' @param link \code{data.frame} with exactly two, named, columns, or a 
#' \code{list} of several such \code{data.frame}s 
#' @returns colnames of `data.frame`(s) in `link`.
#' @noRd
#' 
colnames_link <- function(link){
  if(is.data.frame(link)) return(list(link = names(link)) )
  return(lapply(link, names) )
}