#' Make a MultiFactor
#' @name MultiFactor
#' @rdname MultiFactor
#' @aliases asMultiFactor
#' @description
#' Construct an `MultiFactor` object from a from a validly shaped data frame or
#' list of such data frames.
#' @param x `any` object, most likely `list` of data frames.
#' @param levels an optional named list of vectors of the unique values (as
#'     character strings) that x might have taken. The default is the unique set
#'     of values taken by lapply(x, as.character), sorted into increasing order
#'     of x.
#' @param drop.unmatched `Logical scalar` If `TRUE` (Default), for feature types
#'     that are seen at least twice, exclude features that only present in one
#'     of their respective link data frames.
#' @export
#' @seealso \itemize{
#' \item [kegg_link()]: for an example of valid input.
#' \item [MultiFactor-class()]: for class.
#' \item [MultiFactor-methods()] for methods.
#'}
#' @examples
#' MultiFactor( kegg_link( ) )
#'
MultiFactor <- function(x, levels = NULL, drop.unmatched = TRUE) {
  if(validLinkDF(x)) x <- list(x = x)
  stopifnot( "Input not correctly formatted." = all (
    vapply(as.list(x, use.names = FALSE), validLinkDF, NA, USE.NAMES = FALSE)
    ))

  if(is(x, "MultiFactor")) {
    if(is.null(levels)) {levels <- levels(x)};  x <- x@index
  }

  if(drop.unmatched) { x <- trimMultiFactor(x) }

  m <- mapMultiFactor(x)

  if(is(x, "MultiFactor") && is.null(levels)) {
    levels <- levels(x)
    x <- x@index
    }

  # Integer DF Input
  if(all( vapply(x, validIntLinkDF, NA, USE.NAMES = FALSE))) {
    stopifnot("Input is integers, levels must be provided. " =
                !is.null(levels))
  # Factor DF input
  } else
    if(all( vapply(x, validFactLinkDF, NA, USE.NAMES = FALSE))) {
      if(is.null(levels)) {levels <- factorInputMultiFactorLevels(x, m)}
      x <- listFactRefactor(x, m, levels)
      x <- lapply(x, factToIntDF)
  # Character DF input
    } else
      if(all( vapply(x, validCharLinkDF, NA, USE.NAMES = FALSE))) {
        if(is.null(levels)) {levels <- generateMultiFactorLevels(x, m) }
        x <- listCharToIntegers(x, m, levels)
        }
  # Get rid of row.names.
  x <- lapply(x, `row.names<-.data.frame`, value = NULL)

  # Out
  out <- new("MultiFactor", index = x, levels = levels, map = m)
  validObject(out)
  return(out)
}

#' @noRd
#' @importFrom Matrix sparseMatrix
#' @description
#' Helper function that takes an MultiFactor and returns a sparse biadjacency
#' Matrix with link df names as rownames and id names as colnames. Called
#' internally.
#' @param x a named list of data frames with named character columns.
#' @param int `Logical scalar` Whether to return integer counts (Default) or
#'     otherwise a sparse pattern Matrix.
#' @returns a sparse biadjacency Matrix with link df names as rownames and id
#'     names as colnames. Values count unique features in that position.
#'
mapMultiFactor <- function(x, int = TRUE) {
  # Some flexibility in input
   if( is(x, "MultiFactor") ) {
    x <- x@index }

  all_names <- lapply(x, names)
  i <- factor(rep(names(all_names),
                  vapply(all_names, length, 1, USE.NAMES = FALSE)),
              levels = names(all_names))
  j <- factor(unlist(all_names, use.names = FALSE),
              levels = unique(unlist(all_names, use.names = FALSE)))
  if(!int){
    return( sparseMatrix(i = i, j = j, dimnames = list(levels(i), levels(j))) )
  }
  # Otherwise, add counts.
  mx <- unlist(
    lapply(x, function(y)
      lapply(y, function(z) length(unique(z)))),
    use.names = TRUE)
  return(
    sparseMatrix(i = i, j = j, x = mx, dimnames = list(levels(i), levels(j)))
    )
}


#' @noRd
#' @param x a list in `MultiFactor` formatting.
#' @description Called by `MultiFactor()` if `drop.unmatched` is `TRUE`. Runs
#'     part of `subset` method.
#' @returns a subsetted list with `MultiFactor`.formatting.
#' @importFrom Matrix colSums
#'
trimMultiFactor <- function(x) {
    # Determine positions of feature names that occur in several edge link dfs
    m        <- mapMultiFactor(x, int = FALSE)
    jj       <- colnames(m)[Matrix::colSums(m) > 1]

    # Sequentially subset over feature names
    for(j in jj) {
        # Select all those data frames where that term is mentioned
        ii <- rowsWithCol(m, j, FALSE)
        keep    <- Reduce(intersect, lapply(x[ii], `[[`, j))

        # Filter feature ids in each df to only universally shared ones.
        x[ii] <- lapply(x[ii], function(df) return( df[ df[[j]] %in% keep, ] ))
        }
    return(x)
}


#' @noRd
#' @description
#' helper function to make levels named list for MultiFactor. Not intended to be
#' called directly.
#' @details Especially for large input with many repeated features, it can be
#'     much more efficient to provide levels if known.
#' @param x a named list of data frames with named character columns.
#' @param m Matrix resulting from `mapMultiFactor(x)`
#' @returns a named list of levels.
#'
generateMultiFactorLevels <- function(x, m) {
  lv_names <- colnames(m)
  lev_out   <- lapply(
    seq_along(lv_names),
    function(y) unique(
      unlist(
        lapply(x[rowsWithCol(m, y, FALSE)],
               function(z) unique(z[[ lv_names[y] ]])),
        FALSE, FALSE)
    )
  )
  names(lev_out) <- lv_names
  lev_out
}

#' @noRd
#' @description
#' helper function to make levels named list for MultiFactor. Not intended to be
#' called directly.
#' @details Especially for large input with many repeated features, it can be
#'     much more efficient to provide levels if known.
#' @param x a named list of data frames with named character columns.
#' @param m Matrix resulting from `mapMultiFactor(x)`
#' @returns a named list of levels.
#'
factorInputMultiFactorLevels <- function(x, m) {
  lv_names <- colnames(m)
  lev_out   <- lapply(
    seq_along(lv_names),
    function(y) unique(
      unlist(
        lapply(x[rowsWithCol(m, y, FALSE)],
               function(z) levels(z[[ lv_names[y] ]])),
        FALSE, FALSE)
    )
  )
  names(lev_out) <- lv_names
  lev_out
}


#' @noRd
#' @description
#' Converts input relational data for MultiList into respective from character
#' to integer, according to ranking in overarching levels.
#' @param x list of character data frames
#' @param m map matrix of x
#' @param l levels(x)
#' @returns list of integer data frames accordng to levels.
#'
listCharToIntegers <- function(x, m, l) {

  lv_names <- colnames(m)
  for(id in seq_along(lv_names)) {
    idx <- rownames(m)[m[, id] != 0]

    x[idx] <- lapply(x[idx], charToIntDF, id = lv_names[id], r = l[[id]])
  }
  x
}

#' @noRd
#' @description
#' Converts input relational data for MultiList into respective from character
#' to integer, according to ranking in overarching levels.
#' @param x list of character data frames
#' @param m map matrix of x
#' @param l levels(x)
#' @returns list of integer data frames accordng to levels.
#'
listFactRefactor <- function(x, m, l) {
  lv_names <- colnames(m)
  for(i in seq_along(lv_names)) {
    idx <- rownames(m)[m[, i] != 0]
    x[idx] <- lapply(x[idx], factorToMF, id = lv_names[i], r = l[[i]])
  }
  x
}


#' @param x `data frame`, input
#' @param id `Character scalar`, name of column in x.
#' @param r reference levels
#' @noRd
#'
charToIntDF <- function(x, id, r) {
  x[[ id ]] <- match(x[[ id ]], r)
  return(x)
}

#' @param x `data frame`, input
#' @param id `Character scalar`, name of column in x.
#' @param r reference levels
#' @importFrom forcats lvls_expand
#' @noRd
#'
factorToMF <- function(x, id, r) {
  x[[id]] <- forcats::lvls_expand(x[[id]], r)
  return(x)
}

#' @param x `data frame`, input
#' @noRd
#'
factToIntDF <- function(x) {
  x[] <- lapply(x, as.integer)
  return(x)
}

#' @noRd
#' @description not intended for direct use.
#' @param id feature name, one of `colnames(x)`.
#' @param x named list of data frames with `id %in% colnames()` of those data
#'     frames.
#'
lv_list_char <- function(id, x) sort(
  unique(unlist(lapply(x, function(y)
    unique(y[[id]])), recursive = FALSE, use.names = FALSE))
)


#' @noRd
#' @description
#' Based on base::factor object validation.
#'
validLevels <-  function(levs) {
  if (any(vapply(
    levs, function(x) any(!is.character(x)), NA, USE.NAMES = FALSE
  ))) return("factor levels must be \"character\"")
  if (any(d <- as.logical(vapply(levs, anyDuplicated, 1, USE.NAMES = FALSE))))
    return(cat("duplicated factor levels in level number(s)", which(d)))
  ## 'else'	ok :
  TRUE
}

#' Is this a data.frame with at least two columns, that all are named?
#' @noRd
validLinkDF <- function(x) is.data.frame(x) &&
  NCOL(x) >= 2L && length(colnames(x)) == NCOL(x)

#' @noRd
#' @description Based on `base::factor` object validation.
#' @returns `Logical scalar`, TRUE if valid.
#'
validIntLinkDF <- function(x) validLinkDF(x) &&
  all(vapply(x, is.numeric, NA, USE.NAMES = FALSE))

#' @noRd
#' @description Based on `base::factor` object validation.
#' @returns `Logical scalar`, TRUE if valid.
#'
validCharLinkDF <- function(x) validLinkDF(x) &&
  all(vapply(x, is.character, NA, USE.NAMES = FALSE))

#' @noRd
#' @description Based on `base::factor` object validation.
#' @returns `Logical scalar`, TRUE if valid.
#'
validFactLinkDF <- function(x) validLinkDF(x) &&
  all(vapply(x, is.factor, NA, USE.NAMES = FALSE))

#' @rdname MultiFactor
#' @export
#'
asMultiFactor <- MultiFactor
