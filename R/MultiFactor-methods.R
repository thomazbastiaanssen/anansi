#' S4 Methods for MultiFactor
#' @name MultiFactor-methods
#' @description
#' a set of methods to work with `MultiFactor`, a special type of list.
#'
NULL


#' S4 Methods for MultiFactor
#' @description
#' `getEdgeList`: Return a data frame in edge list format.
#' @rdname MultiFactor-methods
#' @export
#'
setMethod("getEdgeList", "MultiFactor",
          function(x) as.data.frame(do.call(rbind, names(x))))


#' S4 Methods for MultiFactor
#' @description `dim`: Display a vector of dims (n objects, n ids).
#' @rdname MultiFactor-methods
#' @export
#'
setMethod("dim", "MultiFactor", function(x)
  dim(x@map)
)

#' S4 Methods for MultiFactor
#' @description `names`: Display a vector of dims (n objects, n ids).
#' @rdname MultiFactor-methods
#' @export
#'
setMethod("names", "MultiFactor", function(x)
  `names<-`(lapply(x@index, names), rownames(x))
)

#' S4 Methods for MultiFactor
#' @description `dimnames`: Display a vector of dims (n objects, n ids).
#' @rdname MultiFactor-methods
#' @export
#'
setMethod("dimnames", "MultiFactor", function(x)
  dimnames(x@map)
)

#' S4 Methods for MultiFactor
#' @description `[`: Subset based on [rownames(),colnames(x)]
#' @export
#'
setMethod("[", c("MultiFactor", "ANY", "ANY"), definition = function(
    x, i, j, ..., return.list = TRUE, drop = TRUE)
  {
  if ( missing(i) && missing(j) )  return(x)
  d <- dictionary(x)
  x <- unfactor(x)
  if (!missing(i)) ii <- rownames(d[i, , drop = FALSE])
  if (!missing(j)) jj <- colnames(d[, j, drop = FALSE])

  if ( missing(i) ) {
      ii <- rowsWithCol(d, jj, FALSE)
      x  <- lapply(x[ii], `[`, i = jj)

  } else if ( missing(j) ) {
      x <- x[ii]

  } else {
      x  <- lapply(x[ii], `[`, i = jj)
  }

  if(return.list) return(x)

  MultiFactor(x)
})


setReplaceMethod("[", c("MultiFactor", "ANY", "ANY", "list"), def = function(
    x, i, j, ..., value) {
  if (missing(i) && missing(j)) return(value)
  d <- dictionary(x)
  if (!missing(i)) ii <- rownames(d[i ,   , drop = FALSE])
  if (!missing(j)) jj <- colnames(d[  ,  j, drop = FALSE])

  if ( missing(j) ) {
    x@index[ii] <- value
    validObject(x)
    return(x)
  }

  if ( missing(i) ) { ii <- rowsWithCol(d, jj, names = TRUE) }

  for(i in ii) {
      for(j in jj) {
        x@index[[i]][,j] <- value[[i]][,j]
      }
  }
  validObject(x)
  (x)
})




#' @export
#'
setMethod("[[", c("MultiFactor", "ANY"), function(x, i, ...) {
  x@index[[i, ...]]
})

#' @export
#'
setReplaceMethod("[[", c("MultiFactor", "ANY", "ANY"),
                 function(x, i, ..., value) {

                   x@index[[i, ...]] <- value
                   validObject(x)
                   x
                   }
                 )



#' S4 Methods for MultiFactor
#' @description `show`: Display the object
#' @importFrom methods show
#' @importFrom Matrix sparseMatrix printSpMatrix
#' @inheritParams methods::show
#' @rdname MultiFactor-methods
#' @export
#'
setMethod("show",  "MultiFactor", function(object) {
  cat("A list of class ", class(object), ",\n    ",
      NCOL(object), " feature types across ", NROW(object),
      " edge lists.\n\n", sep = "")
  printSpMatrix( object@map )

  cat("\nValues represent unique feature names in that edge list.\n\n",
      "Levels:\n\n", sep = ''
  )
  id_w <- max(nchar(colnames(object)))
  nm_w <- max(nchar(nlevels(object)))
  for(id in colnames(object)) {
    num_lvs <- length(levels(object)[[id]])
    cat(format(id, width = id_w), " : ",
        format(num_lvs, width = nm_w), " Levels: ", sep = "" )

    if(num_lvs > 4L)
      cat(
        levels(object)[[id]][1], levels(object)[[id]][2], "...",
        levels(object)[[id]][num_lvs], "\n",sep = " "
      ) else
        cat(levels(object)[[id]], "\n", sep = " ")

  }
  invisible(NULL)
})

#' @rdname MultiFactor-methods
#' @description Analogous to `factors`. `unfactor(MultiFactor)` returns a named
#'     list with character data frames with same dimensions as input.
#' @importMethodsFrom S4Vectors unfactor
#' @inheritParams S4Vectors::unfactor
#' @returns A named character list
#' @export
#'
setMethod("unfactor", "MultiFactor", function(x) {
  lv <- levels(x)
  ns <- rownames(x)
  x  <- x@index

  x[] <- lapply(x, function(df) {
    for(id in names(df)) {
    df[,id] <- lv[[id]][df[,id]]
    }
    return(df)} )

  return(x)
})

#' @rdname MultiFactor-methods
#' @description Analogous to `factors`. `droplevels(MultiFactor)` returns a
#'     `MultiFactor` with unused levels removed.
#' @importMethodsFrom S4Vectors droplevels
#' @inheritParams base::droplevels
#' @returns A MultiFactor
#' @export
#'
droplevels.MultiFactor <- function(x, ...) {
  lvs   <- levels(x)
  d     <- dictionary(x)
  x.int <- x@index
  for(lv in names(lvs)) {
    rs        <- rowsWithCol(d, lv, names = TRUE)
    x_index   <- lapply(x.int[rs], `[[`, lv )

    x_tot     <- unique(unlist(x_index, use.names = FALSE))

    x_nlevels <- length(lvs[[lv]])
    lvl_ranks <- seq_len(x_nlevels)

    keep_ix   <- which(lvl_ranks %in% x_tot)

    lvs[[lv]] <- lvs[[lv]][keep_ix]
    for(r in rs) { x.int[[r]][,lv]  <- match(x_index[[r]],  table = keep_ix) }
  }
  MultiFactor(x.int, levels = lvs)
}

#' @rdname MultiFactor-methods
#' @export
#'
setMethod("droplevels", "MultiFactor", function(x, ...)
  droplevels.MultiFactor(x, ...))


#mergeROWS()

#' S3/S4 combo for levels.
#' @export
#' @description
#' get object levels
#' @returns a named list of character vectors.
#' @rdname MultiFactor-methods
#'
levels.MultiFactor <- function(x) x@levels

#' @export
#' @rdname MultiFactor-methods
setMethod("levels",  "MultiFactor", levels.MultiFactor)

#' @export
#' @rdname MultiFactor-methods
#' @param value a replacement character vector of suitable dimensions.
#'
setReplaceMethod("levels", "MultiFactor",
                 function(x, value) {
                   x@levels <- value
                   x   } )

#' @export
#' @description
#' get object map
#' @param x `MultiFactor` object
#' @returns a named sparse biadjacency matrix of dimensions (`dimnames(x)`)
#' @rdname MultiFactor-methods
#'
setMethod("dictionary",  "MultiFactor", function(x, ...) x@map)

#' @export
#' @rdname MultiFactor-methods
setReplaceMethod("dictionary", "MultiFactor",
                 function(x, ..., value) {
                   x@map <- value
                   validObject(x)
                   x   } )


#' S4 Methods for MultiFactor
#' @rdname MultiFactor-methods
#' @details
#' `subset`: For `MultiFactor objects`, sub-setting is only applied
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
#' l <- asMultiFactor(kegg_link())
#'
#' # Sub-setting is only performed on data frames that contain the arguments
#' str(subset(x = l, cpd %in% c("C00001", "C00002")))
#'
#' # Several data frames at the same time:
#' subset(x = l, ec %in% c("1.2.3.4", "4.3.2.1"))
#'
setMethod("subset", "MultiFactor", function(x, subset, select, ...) {
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


#' @param d `MultiFactor@map`
#' @param id `Character or Integer scalar`. Selects column of `d`.
#' @param names Whether to return characters (Default) or integer indices.
#' @returns A vector indicating which elements of `MultiFactor` contain `id`.
#' @importFrom Matrix which
#' @noRd
#' @description Helper function for `MultiFactor` to get names or indices of
#' data frames that contain an id column
#'
rowsWithCol <- function(d, id, names = TRUE) {
  rowInds <- Matrix::which((d[,id, drop = FALSE] != 0L))
  if(names){
    rowInds <- rownames(d)[rowInds]
  }
  return(rowInds)
}

#' @param d `MultiFactor@map`
#' @param id `Character or Integer scalar`. Selects row of `d`.
#' @param names Whether to return characters (Default) or integer indices.
#' @returns A vector indicating which feature types are in element `id`.
#' @importFrom Matrix which
#' @noRd
#' @description Helper function for `MultiFactor` to get names or indices of
#'     features contained in a given data frame element of `MultiFactor`.
#'
colsWithRow <- function(d, id, names = TRUE) {
  colInds <- Matrix::which((d[id, , drop = FALSE] != 0L))
  if(names){
    colInds <- colnames(d)[colInds]
  }
  return(colInds)
}
