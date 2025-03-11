#' Accessing and modifying information in AnansiWeb S4 class
#' @name AnansiWeb-methods
#' @description `AnansiWeb` supports `$` operator for getting and
#' assigning values.
#'
#' ` dimnames( x ) ` is shorthand for `dimnames( x$dictionary )` and
#' `names( x )` is in turn shorthand for `names( dimnames(x) )`.
#'
#' @returns a specified `AnansiWeb` object.
#'
#' @seealso \itemize{
#' \item [AnansiWeb-class()].
#' \item [weaveWeb()]: for general use.
#'}
#' @importFrom methods slotNames slot slot<-
#' @examples
#' # prepare an AnansiWeb
#' w <- weaveWeb(cpd ~ ko, link = kegg_link())
#'
#' w$dictionary
#'
#' names(w)
#'
NULL

#' @noRd
#' @export
#' @importFrom utils .DollarNames
.DollarNames.AnansiWeb <- function(x, pattern = "")
  grep(pattern, slotNames(x), value = TRUE)

#' @exportMethod $
#' @inheritParams base::`$`
#' @rdname AnansiWeb-methods
#'
setMethod("$", "AnansiWeb", definition = function(x, name) slot(x, name) )

#' @exportMethod $<-
#' @inheritParams base::`$<-`
#' @rdname AnansiWeb-methods
#'
setReplaceMethod("$", "AnansiWeb", def = function(x, name, value) {
  slot(x, name) <- value
  return(x)}
 )

#' @description `show`: Display the object
#' @importFrom methods show
#' @inheritParams methods::show
#' @rdname AnansiWeb-methods
#' @export
#'
setMethod("show",  "AnansiWeb", function(object) {
  cat(class(object), " object with ", NROW(object$tableX), " observations:\n",
  "    Tables: ", names(object)[1], " (", NROW(object), " features) and ",
  names(object)[2], " (", NCOL(object), " features)\n",
      sep = "")
  cat("Access content with $ operator. ",
      "Collapse with as.list().", sep = "")
  invisible(NULL)
})

#' @rdname AnansiWeb-methods
#' @inheritParams base::dimnames
#' @export
#'
setMethod("dimnames", "AnansiWeb",
          function(x) dimnames(x@dictionary)
)

#' @rdname AnansiWeb-methods
#' @inheritParams base::dim
#' @export
#'
setMethod("dim", "AnansiWeb",
          function(x) dim(x@dictionary)
)

#' @rdname AnansiWeb-methods
#' @inheritParams base::names
#' @export
#'
setMethod("names", "AnansiWeb", function(x) names( dimnames( x@dictionary) ))

#' @noRd
#'
tell_F <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@f.values)
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_F))
  }
}

#' @noRd
#'
tell_T <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@t.values)
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_T))
  }
}

#' @noRd
#'
tell_P <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@p.values)
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_P))
  }
}


#' @noRd
#'
tell_e <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@estimates)
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_e))
  }
}


#' @noRd
#'
tell_df1 <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@df[1])
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_df1))
  }
}

#' @noRd
#'
tell_df2 <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@df[2])
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_df2))
  }
}

#' @noRd
#'
tell_dfr <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@df[3])
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_dfr))
  }
}

################################################################################
################################################################################

#' Is this a data.frame with exactly two columns that are named?
#' @noRd
validWeb <- function(x) {
  y_names  <- identical(rownames(x), colnames(x$tableY))
  x_names  <- identical(colnames(x), colnames(x$tableX))
  s_names  <- identical(rownames(x$tableY), rownames(x$tableX))
  meta_dim <- any(NROW(x$metadata) == NROW(x$tableY),
                  prod(dim(x$metadata)) <= 1)
  if(!y_names)
    message("colnames(tableY), rownames(dictionary) not identical.")
  if(!x_names)
    message("colnames(tableX), colnames(dictionary) not identical.")
  if(!s_names)
    message("rownames(tableX), rownames(tableX) not identical.")
  if(!meta_dim)
    message("NROW(metadata) does not equal rows of tableY, tableX.")

  return(all(y_names, y_names, s_names, meta_dim))
}
