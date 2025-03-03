#' Accessing and modifying information in anansiWeb S4 class
#' @name anansiWeb-methods
#' @description \code{anansiWeb} supports \code{$} operator for getting and 
#' assigning values.
#' 
#' \code{ dimnames( x ) } is shorthand for \code{dimnames( x$dictionary )} and 
#' \code{terms( x )} is in turn shorthand for \code{names( dimnames(x) )}. 
#' 
#' @returns a specified \code{anansiWeb} object. 
#' 
#' @seealso \itemize{
#' \item \code{\link{anansiWeb-class}}. 
#' \item \code{\link{weaveWeb}}: for general use.
#'}
#' @importFrom methods slotNames slot slot<- 
#' @examples
#' # prepare an anansiWeb
#' w <- weaveWeb(cpd ~ ko)
#' 
#' w$dictionary
#' 
#' terms(w) 
#'   
NULL

#' @noRd
#' @export
#' @importFrom utils .DollarNames
.DollarNames.anansiWeb <- function(x, pattern = "")
  grep(pattern, slotNames(x), value = TRUE)

#' @exportMethod $
#' @inheritParams base::`$`
#' @rdname anansiWeb-methods
#' 
setMethod("$", "anansiWeb", definition = function(x, name) slot(x, name) )

#' @exportMethod $<-
#' @inheritParams base::`$<-`
#' @rdname anansiWeb-methods
#' 
setReplaceMethod("$", "anansiWeb", def = function(x, name, value) {
  slot(x, name) <- value
  return(x)}
 )

#' @rdname anansiWeb-methods
#' @export
#' 
setMethod("dimnames", "anansiWeb", 
          function(x) dimnames(x@dictionary)
)

#' @rdname anansiWeb-methods
#' @export
#' 
setMethod("names", "anansiWeb", 
          function(x) names( dimnames( x@dictionary) ) 
)

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