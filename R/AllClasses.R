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
#' @importClassesFrom S4Vectors DataFrame
#' @seealso \itemize{
#' \item [AnansiWeb()]: for general use.
#' \item [AnansiWeb-methods()] for methods, including `$`
#' operator.
#'}
#'
setClass("AnansiWeb",
         slots = c(
           tableY     = "matrix",
           tableX     = "matrix",
           dictionary = "Matrix",
           metadata   = "DataFrame"
         )
)

#' AnansiLinkMap S4 container class
#' @name AnansiLinkMap-class
#' @description
#' `AnansiLinkMap` is an S4 class containing one or several data frames
#' structured as edge lists from the `igraph` package.
#' @export
#' @seealso \itemize{
#' \item [AnansiLinkMap()]: for general use.
#' \item [AnansiLinkMap-methods()] for methods
#' \item [igraph::igraph()].
#'}
#'
setClass("AnansiLinkMap",
         contains = "list")

#' is valid AnansiLinkMap?
#' @noRd
#' @description
#' returns TRUE if input is in the right format to be an AnansiLinkMap object
#' @param object
#' `any` object, but not much will happen unless the object's class has a
#' formal definition.
#' @importFrom methods validObject
#' @returns `TRUE` if passes, character vector otherwise.
#'
setValidity("AnansiLinkMap", method = function(object) ifelse(
    test = validLinkDF(object) || all(unlist(lapply(object, validLinkDF))),
    yes = TRUE,
    no = "object is not in a valid format.")
)

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
