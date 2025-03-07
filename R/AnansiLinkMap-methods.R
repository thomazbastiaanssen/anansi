#' S4 Methods for AnansiLinkMap
#' @name AnansiLinkMap-methods
#' @description
#' a set of methods to work with `anansiLinkWeb`, a special type of list.
#'
NULL


#' S4 Methods for AnansiLinkMap
#' @description
#' `getEdgeList`: Return a data frame in edge list format.
#' @rdname AnansiLinkMap-methods
#' @export
#'
setMethod("getEdgeList", "AnansiLinkMap",
          function(x) as.data.frame(do.call(rbind, names(x))))

#' S3/S4 combo for levels.AnansiLinkMap
#' @description
#' following [S4Vectors:levels.Rle()]
#' @rdname AnansiLinkMap
#'
levels.AnansiLinkMap <- function(x) lapply(colnames(x), lv_list_char, x)
setMethod("levels", "AnansiLinkMap", levels.AnansiLinkMap)

setReplaceMethod("levels", "AnansiLinkMap",
                 function(x, value) {
                   levels(x@values) <- value
                   validObject(x)
                   x   } )
#' S4 Methods for AnansiLink
#' Map
#' @description `names`: Display a list of column names from AnansiLinkMap
#' @rdname AnansiLinkMap-methods
#' @export
#'
setMethod("names", "AnansiLinkMap", function(x) lapply(x, names) )

#' S4 Methods for AnansiLinkMap
#' @description `dimnames`: Display a vector of id names.
#' @rdname AnansiLinkMap-methods
#' @export
#'
setMethod("dimnames", "AnansiLinkMap", function(x)
  list(
    names(names(x)),
    unique(unlist(names(x), use.names = FALSE))
  )
)

#' S4 Methods for AnansiLinkMap
#' @description `dim`: Display a vector of dims (n objects, n ids).
#' @rdname AnansiLinkMap-methods
#' @export
#'
setMethod("dim", "AnansiLinkMap", function(x)
  c(length(x),length(unique(unlist(names(x), use.names = FALSE))))
  )

#' S4 Methods for AnansiLinkMap
#' @description `show`: Display the object
#' @importFrom methods show
#' @importFrom Matrix sparseMatrix printSpMatrix
#' @inheritParams methods::show
#' @rdname AnansiLinkMap-methods
#' @export
#'
setMethod("show",  "AnansiLinkMap", function(object) {
    cat("A list of class ", class(object), ",\n    ",
        NCOL(object), " feature types across ", NROW(object), " edge lists.\n\n", sep = "")
    c.n <- colnames(object)
    lv_lst <- lapply(c.n, lv_list_char, x = object)

  i <- factor(rep(rownames(object), each = 2))
  j <- factor(unlist(names(object), use.names = FALSE))
  x <- unlist(
    lapply(object, function(x)
      lapply(x, function(y)
        length(unique(y))
        )), use.names = FALSE)
 sm <- sparseMatrix(
   i = i,
   j = j,
   x = x,
   dimnames = dimnames(object))
 printSpMatrix(sm)

 cat("\nValues represent unique feature names in that edge list.")
 invisible(NULL)
})


#' S4 Methods for AnansiLinkMap
#' @rdname AnansiLinkMap-methods
#' @details
#' `subset`: For `AnansiLinkMap objects`, sub-setting is only applied
#' to data frames compatible with the expression. The rest are returned
#' unaltered. Modeled after `subset()`.
#'
#' @param subset
#' `logical expression` indicating rows to keep. Must contain variables
#' found as column names.
#' @param select `expression`. Which column names to consider. If missing
#' (Default), consider all column names.
#' @inheritParams BiocGenerics::subset
#' @importMethodsFrom BiocGenerics subset
#' @export
#' @seealso [BiocGenerics::subset()].
#' [weaveWeb()] for the AnansiWeb constructor functions that
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
setMethod("subset", "AnansiLinkMap", function(x, subset, select, ...) {
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

#' @noRd
#' @description not intended for direct use.
#' @param id feature name, one of `colnames(x)`.
#' @param x `AnansiLinkWeb`
#'
lv_list_char <- function(id, x)
  sort(unique(unlist(lapply(x, function(y)
    unique(as.character(y[[id]]))), recursive = FALSE, use.names = FALSE)))

#' @noRd
#' @description not intended for direct use.
#' @param id feature name, one of `colnames(x)`.
#' @param x `AnansiLinkWeb`
#'
lv_list_factor <- function(id, x)
  unique(unlist(lapply(x, function(y)
    unique(levels(y[[id]]))), recursive = FALSE, use.names = FALSE))

#' Is this a data.frame with exactly two columns that are named?
#' @noRd
validLinkDF <- function(x) is.data.frame(x) &&
  NCOL(x) == 2L &&
  length(colnames(x)) == 2L
