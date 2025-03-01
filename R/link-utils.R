#' make an igraph mapping across the content of link
#' @details Make an \code{igraph} graph from the  content of \code{link}.  
#' @param link 
#' \code{data.frame} with exactly two, named, columns, or a \code{list} of 
#' several such \code{data.frame}s 
#' @returns an \code{igraph} graph. 
#' @importFrom igraph graph_from_edgelist
#' @seealso \code{\link{weaveWeb}} for the anansiWeb constructor functions that 
#' take link data frames, \code{\link[igraph:igraph]{igraph}} for output object.
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
#' @inheritParams graph_from_link
#' @returns colnames of `data.frame`(s) in `link`.
#' @noRd
#' 
colnames_link <- function(link){
  if(is.data.frame(link)) return(list(link = names(link)) )
  return(lapply(link, names) )
}

#' Subset a valid link object 
#' @details
#' Convenience function for valid link objects. For link lists, sub-setting is 
#' only applied to data frames compatible with the expression. The rest are 
#' returned unaltered. Modeled after \code{subset()}. 
#' 
#' @inheritParams graph_from_link
#' @param subset 
#' \code{logical expression} indicating rows to keep. Must contain variables 
#' found as column names.
#' @export
#' @seealso \code{\link[base:subset]{subset}}, where this function is based on, 
#' \code{\link{weaveWeb}} for the anansiWeb constructor functions that 
#' take link data frames.
#' @examples
#' # Works for single link data frame, as well as list of multiple such objects.
#' head(subset_link(link = ec2cpd, cpd %in% c("C00001", "C00002")))
#' 
#' # Sub-setting is only performed on data frames that contain the arguments
#' str(subset_link(link = kegg_link(), cpd %in% c("C00001", "C00002")))
#' 
#' # Several data frames at the same time:
#' str(subset_link(link = kegg_link(), ec %in% c("1.2.3.4", "4.3.2.1")))
#' 
subset_link <- function(link, subset){
  
  stopifnot("'link' must be one or several data.frame objects." =
              validLink(link))
  
  # Mimic missing behaviour of base::subset 
  if(missing(subset)) subset <- quote(TRUE) else subset <- substitute(subset)
  
  # Simple data frame case before anything else. 
  if(validLinkDF(link)) { 
    r <- eval( subset, link, parent.frame() )
    return(link[r,]) }
  
  # Now handle list, start with selection
  sub.vars <- all.vars(subset)
  l.names <- colnames_link(link)
  
  # Select those data frames where all terms are mentioned
  sub.ind <- unlist(lapply(l.names, function(x) all( sub.vars %in% x )))
  
 # return( )
 link[sub.ind] <- lapply(link[sub.ind], function(x) {
    r <- eval(subset, x, parent.frame() )
    return(x[r,]) })
 
 return(link)

}
