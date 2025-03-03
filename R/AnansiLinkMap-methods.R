#' Make and anansiLinkMap
#' @name LinkMap
#' @rdname LinkMap
#' @aliases anansiLinkMap asLinkMap
#' @description
#' \code{asLinkMap()} constructs an \code{anansiLinkMap} object from a validly 
#' formed data frame or list of such data.frames.  
#' @param x \code{any} object, most likely \code{list} of data frames. 
#' @export
#' @seealso \itemize{
#' \item \code{\link{kegg_link}}: for an example of valid input. 
#' \item \code{\link{anansiLinkMap-class}}: for class.
#' \item \code{\link{anansiLinkMap-methods}} for methods.
#'}
#' @examples
#' asLinkMap( kegg_link() )
#' 
asLinkMap <- function(x) {
  if(validLinkDF(x)) x <- list(link = x)
  
  linkMap <- new("anansiLinkMap", x)
  
  linkMap@edgelist <- as.data.frame(do.call(rbind, names(linkMap)))
  
  return(linkMap)
}

#' Methods for anansiLinkMap S4 class
#' @rdname anansiLinkMap-methods
#' @description
#' \code{clean}: Remove all unpaired feature names from an \code{anansiLinkMap}. 
#' @inheritParams ShortRead::clean
#' @export
#' @returns 
#' \code{clean}: An instance of \code{class(object)}, containing only feature 
#' names that are present in all constituent data frames with that column. 
#' 
setMethod("clean", "anansiLinkMap", 
          function(object, ...) {
            return(  "hue"    )
          })

#' @description \code{names}: Display a list of column names from anansiLinkMap 
#' @rdname anansiLinkMap-methods
#' @export
#'
setMethod("names", "anansiLinkMap", function(x) lapply(x, names) )

#' @rdname anansiLinkMap-methods
#' @details
#' \code{subset}: For \code{anansiLinkMap objects}, sub-setting is only applied 
#' to data frames compatible with the expression. The rest are returned 
#' unaltered. Modeled after \code{subset()}. 
#' 
#' @param subset 
#' \code{logical expression} indicating rows to keep. Must contain variables 
#' found as column names.
#' @inheritParams BiocGenerics::subset
#' @export
#' @seealso \code{\link[BiocGenerics:subset]{subset}}. 
#' \code{\link{weaveWeb}} for the anansiWeb constructor functions that 
#' take link data frames.
#' @examples
#' # prep input
#' l <- asLinkMap(kegg_link())
#' 
#' # Sub-setting is only performed on data frames that contain the arguments
#' str(subset(x = l, cpd %in% c("C00001", "C00002")))
#' 
#' # Several data frames at the same time:
#' str(subset(x = l, ec %in% c("1.2.3.4", "4.3.2.1")))
#' 
setMethod("subset", "anansiLinkMap", 
          function(x, subset, ...) {
            stopifnot("'link' must be one or several data.frame objects." =
                        validLinkMap(x))
            
            # Mimic missing behaviour of base::subset 
            if(missing(subset)) subset <- quote(TRUE) else subset <- substitute(subset)
            
            # Now handle list, start with selection
            sub.vars <- all.vars(subset)
            x.names <- names(x)
            
            # Select those data frames where all terms are mentioned
            sub.ind <- unlist(lapply(x.names, function(y) all( sub.vars %in% y )))
            
            # Subset
            x[sub.ind] <- lapply(x[sub.ind], function(y) {
              r <- eval(subset, y, parent.frame() )
              return(y[r,]) })
            
            return(x)
          }
)

#' Methods for anansiLinkMap S4 class
#' @rdname anansiLinkMap-methods
#' @description
#' \code{clean}: Remove all unpaired feature names from an \code{anansiLinkMap}. 
#' @inheritParams ShortRead::clean
#' @export
#' @returns 
#' \code{clean}: An instance of \code{class(object)}, containing only feature 
#' names that are present in all constituent data frames with that column. 
#' 
setMethod("clean", "anansiLinkMap", 
          function(object, ...) {
            return(  "hue"    )
          })

#' @description \code{names}: Display a list of column names from anansiLinkMap 
#' @rdname anansiLinkMap-methods
#' @export
#'
setMethod("names", "anansiLinkMap", function(x) lapply(x, names) )

#' @rdname anansiLinkMap-methods
#' @details
#' \code{subset}: For \code{anansiLinkMap objects}, sub-setting is only applied 
#' to data frames compatible with the expression. The rest are returned 
#' unaltered. Modeled after \code{subset()}. 
#' 
#' @param subset 
#' \code{logical expression} indicating rows to keep. Must contain variables 
#' found as column names.
#' @inheritParams BiocGenerics::subset
#' @export
#' @seealso \code{\link[BiocGenerics:subset]{subset}}. 
#' \code{\link{weaveWeb}} for the anansiWeb constructor functions that 
#' take link data frames.
#' @examples
#' # prep input
#' l <- asLinkMap(kegg_link())
#' 
#' # Sub-setting is only performed on data frames that contain the arguments
#' str(subset(x = l, cpd %in% c("C00001", "C00002")))
#' 
#' # Several data frames at the same time:
#' str(subset(x = l, ec %in% c("1.2.3.4", "4.3.2.1")))
#' 
setMethod("subset", "anansiLinkMap", 
          function(x, subset, ...) {
            stopifnot("'link' must be one or several data.frame objects." =
                        validLinkMap(x))
            
            # Mimic missing behaviour of base::subset 
            if(missing(subset)) subset <- quote(TRUE) else subset <- substitute(subset)
            
            # Now handle list, start with selection
            sub.vars <- all.vars(subset)
            x.names <- names(x)
            
            # Select those data frames where all terms are mentioned
            sub.ind <- unlist(lapply(x.names, function(y) all( sub.vars %in% y )))
            
            # Subset
            x[sub.ind] <- lapply(x[sub.ind], function(y) {
              r <- eval(subset, y, parent.frame() )
              return(y[r,]) })
            
            return(x)
          }
)

#' Get a graph object out of an anansiLinkMap. 
#' @rdname anansiLinkMap-methods
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

setMethod("getGraph", "list", 
          function(x, format = "igraph", ...) getGraph( asLinkMap(x), format)
)

#' Is this a data.frame with exactly two columns that are named?
#' @noRd
validLinkDF <- function(x) is.data.frame(x) && 
  NCOL(x) == 2L &&
  length(colnames(x)) == 2L
