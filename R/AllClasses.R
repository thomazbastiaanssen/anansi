#' anansiWeb S4 container class
#' @name anansiWeb-class
#' @description
#' \code{anansiWeb} is an S4 class containing two feature tables as well as a
#' dictionary to link them. \code{anansiWeb} is the main container that will
#' hold your input data throughout the \code{anansi} pipeline.
#' @slot tableY \code{matrix} of metabolomics data. Rows are samples and columns
#'     are features.
#' @slot tableX \code{matrix} of functional data. Rows are samples and columns
#'     are features.
#' @slot dictionary \code{Matrix}, binary adjacency matrix. Optionally sparse.
#'     Typically generated using the\code{weaveWeb()} function.
#' @importClassesFrom Matrix Matrix
#' @seealso \itemize{
#' \item \code{\link{weaveWeb}}: for general use.
#' \item \code{\link{anansiWeb-methods}} for methods, including \code{$}
#' operator.
#'}
#'
setClass("anansiWeb",
         slots = c(
           tableY     = "matrix",
           tableX     = "matrix",
           dictionary = "Matrix"
         )
)

#' anansiLinkMap S4 container class
#' @name anansiLinkMap-class
#' @description
#' \code{anansiLinkMap} is an S4 class containing one or several data frames
#' structured as edge lists from the \code{igraph} package.
#' @export
#' @seealso \itemize{
#' \item \code{\link{LinkMap}}: for general use.
#' \item \code{\link{anansiLinkMap-methods}} for methods
#' \item \code{\link[igraph:igraph]{igraph}}.
#'}
#'
setClass("anansiLinkMap",
         contains = "list")

#' is valid anansiLinkMap?
#' @noRd
#' @description
#' returns TRUE if input is in the right format to be an anansiLinkMap object
#' @param object
#' \code{any} object, but not much will happen unless the object's class has a
#' formal definition.
#' @importFrom methods validObject
#' @returns \code{TRUE} if passes, character vector otherwise.
#'
setValidity("anansiLinkMap", method = function(object) ifelse(
    test = validLinkDF(object) || all(unlist(lapply(object, validLinkDF))),
    yes = TRUE,
    no = "object is not in a valid format.")
)

#' An S4 class to contain all \code{anansi} stats results so that they can
#' easily be extracted.
#'
#' @slot subject A character that describes the data that was queried.
#' @slot type A character that describes type of parameter contained in the
#'     \code{estimates} slot. For example r.values for correlations or r.squared
#'     for models.
#' @slot df a vector of length 2, containing df1 and df2 corresponding to the
#'     F-ratio considered.
#' @slot estimates A matrix containing the estimates for the parameters named in
#'     the \code{type} slot.
#' @slot f.values A matrix containing the f-values, for least-squares.
#' @slot t.values A matrix containing the t-values, for correlations.
#' @slot p.values A matrix containing the p.values for the parameters named in
#'     the \code{type} slot.
#' @description \code{anansiTale} is the main container that will hold your
#'     stats output data coming out of the \code{anansi} pipeline.
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
