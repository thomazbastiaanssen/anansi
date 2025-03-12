#' Make aa MultiFactor
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
#' @export
#' @seealso \itemize{
#' \item [kegg_link()]: for an example of valid input.
#' \item [MultiFactor-class()]: for class.
#' \item [MultiFactor-methods()] for methods.
#'}
#' @examples
#' MultiFactor( kegg_link( ) )
#'
MultiFactor <- function(x, levels) {
  if(validLinkDF(x)) x <- list(x = x)
  stopifnot( "Input not correctly formatted." =
               all( vapply(x, validLinkDF, NA, USE.NAMES = FALSE)))

  m <- mapMultiFactor(x)

  if(all( vapply(x, validIntLinkDF, NA, USE.NAMES = FALSE))) {
    stopifnot(
      "Input is integers, levels must be provided. " = !missing(levels)
    )
  } else
    if(all( vapply(x, validFactLinkDF, NA, USE.NAMES = FALSE))) {
      if(missing(levels)) levels <- factorInputMultiFactorLevels(x, m)
      x <- listFactRefactor(x, m, levels)
      x <- lapply(x, factToIntDF)
    } else
      if(all( vapply(x, validCharLinkDF, NA, USE.NAMES = FALSE))) {
        if(missing(levels)) levels <- generateMultiFactorLevels(x, m)
        x <- listCharToIntegers(x, m, levels)
      }

  out <- new("MultiFactor", x, levels = levels, map = m)

  validObject(out)

  return(out)
}

#' @noRd
#' @importFrom Matrix sparseMatrix
#' @description
#' helper function to make mapping matrix, internal use.
#' @param x a named list of data frames with named character columns.
#' @returns a sparse matrix indicating which ids can be found in which data
#'     frames.
#'
mapMultiFactor <- function(x) {
  all_names <- lapply(x, names)
  i <- factor(rep(names(all_names),
                  vapply(all_names, length, 1, USE.NAMES = FALSE)),
              levels = names(all_names))
  j <- factor(unlist(all_names, use.names = FALSE),
              levels = unique(unlist(all_names, use.names = FALSE)))
  mx <- unlist(
    lapply(x, function(y)
      lapply(y, function(z) length(unique(z)))),
    use.names = TRUE)
  sparseMatrix(
    i = i, j = j, x = mx, dimnames = list(levels(i), levels(j))
  )
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
        lapply(x[names(which(m[, y] != 0 ))],
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
        lapply(x[names(which(m[, y] != 0 ))],
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
  x[[id]] <- lvls_expand(x[[id]], r)
  return(x)
}

#' @param x `data frame`, input
#' @noRd
#'
factToIntDF <- function(x) {
  x <- as.data.frame.list(lapply(x, as.integer))
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
validLevels <-  function(x) {
  levs <- levels(x)
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

#' @rdname MultiFactor
#' @description
#' Helper function that takes an MultiFactor and returns a sparse biadjacency
#' Matrix with link df names as rownames and id names as colnames. Called
#' internally.
#' @returns a sparse biadjacency Matrix with link df names as rownames and id
#'     names as colnames
#' @importFrom Matrix sparseMatrix
#'
linkMatrix <- function(x){
  i <- factor(rep(rownames(x), each = 2))
  j <- factor(unlist(names(x), use.names = FALSE))
  sparseMatrix(i, j, dimnames = dimnames(x))
}


