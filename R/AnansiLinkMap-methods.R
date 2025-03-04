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
  validObject(linkMap)

  return(linkMap)
}

#' S4 Methods for anansiLinkMap
#' @name anansiLinkMap-methods
#' @description
#' a set of methods to work with \code{anansiLinkWeb}, a special type of list.
#'
NULL

#' S4 Methods for anansiLinkMap
#' @description \code{names}: Display a list of column names from anansiLinkMap
#' @rdname anansiLinkMap-methods
#' @export
#'
setMethod("names", "anansiLinkMap", function(x) lapply(x, names) )

#' S4 Methods for anansiLinkMap
#' @description \code{show}: Display the object
#' @importFrom methods show
#' @inheritParams methods::show
#' @rdname anansiLinkMap-methods
#' @export
#'
setMethod("show",  "anansiLinkMap", function(object) {
    cat("An object of class ", class(object), ":\n", sep = "")
    show(lapply(object,
                function(x) (apply(x, 2, function(y) length(unique(y)))))
         )
    invisible(NULL)
})

#' S4 Methods for anansiLinkMap
#' @description
#' \code{getEdgeList}: Return a data frame in edge list format.
#' @rdname anansiLinkMap-methods
#' @export
#'
setMethod("getEdgeList", "anansiLinkMap",
          function(x) as.data.frame(do.call(rbind, names(x))))


#' S4 Methods for anansiLinkMap
#' @rdname anansiLinkMap-methods
#' @details
#' \code{subset}: For \code{anansiLinkMap objects}, sub-setting is only applied
#' to data frames compatible with the expression. The rest are returned
#' unaltered. Modeled after \code{subset()}.
#'
#' @param subset
#' \code{logical expression} indicating rows to keep. Must contain variables
#' found as column names.
#' @param select \code{expression}. Which column names to consider. If missing
#' (Default), consider all column names.
#' @inheritParams BiocGenerics::subset
#' @importMethodsFrom BiocGenerics subset
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
#' subset(x = l, ec %in% c("1.2.3.4", "4.3.2.1"))
#'
setMethod("subset", "anansiLinkMap", function(x, subset, select, ...) {
  validObject(x); x.names <- names(x)
  # PART I: SUBSETTING
  if(!missing(subset)) {
    subset <- substitute(subset); sub.vars <- all.vars(subset)
    # Select those data frames where all terms are mentioned
    sub.ind <- unlist(lapply(x.names, function(y) all( sub.vars %in% y )))
    # Subset them
    x[sub.ind] <- lapply(x[sub.ind], function(y) {
      r <- eval(subset, y, parent.frame() )
      return(y[r,]) })  }
  # Return now if only one df.
  if(length(x) == 1L) return(x)

  # PART II: SELECTING
  if(missing(select)) {
    id.vec   <- unlist(x.names, use.names = FALSE)
    id.share <- id.vec[duplicated(id.vec)]
    sel.vars <- id.share} else sel.vars <- all.vars(substitute(select))
  for(v in sel.vars) {
    # Select those data frames where all terms are mentioned
    s.ind <- unlist(lapply(x.names, function(y) v %in% y ))
    sel.obj <- x[s.ind]
    keep    <- Reduce(intersect, lapply(sel.obj, function(df) df[,v]))
    # Filter feature ids in each df to only include universally shared ones.
    x[s.ind] <- lapply(sel.obj, function(df) return( df[df[[v]] %in% keep,] ))
  }
  return(x)
})

################################################################################
################################################################################

#' Is this a data.frame with exactly two columns that are named?
#' @noRd
validLinkDF <- function(x) is.data.frame(x) &&
  NCOL(x) == 2L &&
  length(colnames(x)) == 2L
