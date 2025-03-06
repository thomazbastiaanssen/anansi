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

#' @rdname AnansiWeb-methods
#' @inheritParams base::dimnames
#' @export
#'
setMethod("dimnames", "AnansiWeb",
          function(x) dimnames(x@dictionary)
)

#' @rdname AnansiWeb-methods
#' @inheritParams base::names
#' @export
#'
setMethod("names", "AnansiWeb", function(x) names( dimnames( x@dictionary) ))

###############################################################################
###############################################################################
#
# setAs

#' Accessing and modifying information in AnansiWeb S4 class
#' @name as
#' @rdname AnansiWeb-methods
#' @importFrom methods as
#' @export
#'
setAs(from = "AnansiWeb", to = "list", def = function(from) {
  out <- list(tableY = from@tableY, tableX = from@tableX,
              dictionary = from@dictionary, metadata = from@metadata)
  names(out)[c(1L, 2L)] <- names(from)
  out
})

#' Accessing and modifying information in AnansiWeb S4 class
#' @name as
#' @rdname AnansiWeb-methods
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#'
setAs(from = "AnansiWeb", to = "MultiAssayExperiment", def = function(from) {

  tY  <- t(from@tableY)
  tX  <- t(from@tableX)
  to_exp <- ExperimentList(
    y = SummarizedExperiment(tY),
    x = SummarizedExperiment(tX)
  )
  names(to_exp) <- names(from)

  to_md  <- list(dictionary = from@dictionary)
  to_cd  <- from@metadata

  MultiAssayExperiment(
    experiments = to_exp,
    metadata = to_md,
    colData = to_cd
    )
})


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
