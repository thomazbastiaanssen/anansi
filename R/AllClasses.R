#' AnansiWeb S4 container class
#' @name AnansiWeb-class
#' @description
#' `AnansiWeb` is an S4 class containing two feature tables as well as a
#' dictionary to link them. `AnansiWeb` is the main container that will
#' hold your input data throughout the `anansi` pipeline.
#' @slot tableY `matrix` of metabolomics data. Rows are samples and columns
#'     are features.
#' @slot tableX `matrix` of functional data. Rows are samples and columns
#'     are features.
#' @slot dictionary `Matrix`, binary adjacency matrix. Optionally sparse.
#'     Typically generated using the`weaveWeb()` function.
#' @importClassesFrom Matrix Matrix
#' @importClassesFrom S4Vectors Annotated
#' @seealso \itemize{
#' \item [AnansiWeb()]: for general use.
#' \item [AnansiWeb-methods()] for methods, including `$`
#' operator.
#'}
#'
setClass("AnansiWeb",
         contains = "Annotated",
         slots = c(
           tableY     = "matrix",
           tableX     = "matrix",
           dictionary = "Matrix",
           metadata   = "list"
         )
)

#' is valid AnansiWeb?
#' @noRd
#' @description
#' returns TRUE if input is in the right format to be an AnansiWeb object
#' @param object
#' `any` object, but not much will happen unless the object's class has a
#' formal definition.
#' @importFrom methods validObject
#' @returns `TRUE` if passes, character vector otherwise.
#'
setValidity("AnansiWeb", method = function(object) ifelse(
  test = validWeb(object),
  yes = TRUE,
  no = "object is not in a valid format.")
)

#' MultiFactor S4 container class
#' @name MultiFactor-class
#' @description
#' `MultiFactor` is an S4 class containing one or several data frames
#' structured as edge lists from the `igraph` package.
#' @export
#' @seealso \itemize{
#' \item [MultiFactor()]: for general use.
#' \item [MultiFactor-methods()] for methods
#' \item [igraph::igraph()].
#'}
#'

#' MultiFactor S4 container class
#' @description
#' `MultiFactor` is an S4 class to manage multiple sets of factors. Methods for
#' `MultiFactor` aim to follow `factor` behaviour.
#' @slot levels `Named list of character vectors`
#' @slot map `(sparse)Matrix` specifying which elements contain which levels.
#' @importClassesFrom Matrix Matrix
#'
#' @export
setClass("MultiFactor",
         contains = "list",
         slots = c(levels  = "list",
                   map     = "Matrix")
)

#' is valid MultiFactor?
#' @noRd
#' @description
#' returns TRUE if input is in the right format to be an MultiFactor object
#' @param object
#' `any` object, but not much will happen unless the object's class has a
#' formal definition.
#' @importFrom methods validObject
#' @returns `TRUE` if passes, character vector otherwise.
#'
setValidity("MultiFactor", method = function(object) ifelse(
  test = validMultiFactor(object),
  yes = TRUE,
  no = "object is not in a valid format.")
)


#' Is this a data.frame with exactly two columns that are named?
#' @noRd
validMultiFactor <- function(x) {

  levels_valid <- validLevels(x)
  values_valid <- vapply(x, validIntLinkDF, NA, USE.NAMES = FALSE)
  no_missing   <- ! any(vapply(x, anyNA, NA, USE.NAMES = FALSE))

  if(!isTRUE(levels_valid))
    message("Levels are not structured correctly. ")
  if(!isTRUE(all(values_valid)))
    message("List content in positions ",
            paste(which(!isTRUE(values_valid)), collapse = ", "),
            " not structured correctly. ")
  if(!no_missing)
    message("Missing values are not allowed.")

  if(
    all(levels_valid,
        isTRUE(values_valid),
        isTRUE(all(no_missing))
    )
  ) return( TRUE )

}

#' An S4 class to contain all `anansi` stats results so that they can
#' easily be extracted.
#'
#' @slot subject A character that describes the data that was queried.
#' @slot type A character that describes type of parameter contained in the
#'     `estimates` slot. For example r.values for correlations or r.squared
#'     for models.
#' @slot df a vector of length 2, containing df1 and df2 corresponding to the
#'     F-ratio considered.
#' @slot estimates A matrix containing the estimates for the parameters named in
#'     the `type` slot.
#' @slot f.values A matrix containing the f-values, for least-squares.
#' @slot t.values A matrix containing the t-values, for correlations.
#' @slot p.values A matrix containing the p.values for the parameters named in
#'     the `type` slot.
#' @description `anansiTale` is the main container that will hold your
#'     stats output data coming out of the `anansi` pipeline.
#'
setClass("anansiTale",
  slots = c(
    subject   = "character",
    type      = "character",
    df        = "numeric",
    estimates = "matrix",
    f.values  = "matrix",
    t.values  = "matrix",
    p.values  = "matrix"
  )
)
