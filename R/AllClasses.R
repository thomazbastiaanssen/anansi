#' anansiWeb S4 container class
#' @name anansiWeb-class
#' @description 
#' \code{anansiWeb} is an S4 class containing two feature tables as well as a 
#' dictionary to link them. \code{anansiWeb} is the main container that will 
#' hold your input data throughout the \code{anansi} pipeline.
#' @slot tableY \code{matrix} of metabolomics data. Rows are samples and columns
#' are features.
#' @slot tableX \code{matrix} of functional data. Rows are samples and columns 
#' are features.
#' @slot dictionary \code{Matrix}, binary adjacency matrix. Optionally sparse. 
#' Typically generated using the\code{weaveWeb()} function.
#' @importClassesFrom Matrix Matrix
#' @seealso \itemize{
#' \item \code{\link{weaveWeb}}: for general use.
#' \item \code{\link{anansiWeb-methods}} for methods, including \code{$} 
#' operator.
#'}
setClass("anansiWeb",
         slots = c(
           tableY     = "matrix",
           tableX     = "matrix",
           dictionary = "Matrix"
         )
)
 
#' An S4 class to contain all \code{anansi} stats results so that they can
#' easily be extracted.
#'
#' @slot subject A character that describes the data that was queried.
#' @slot type A character that describes type of parameter contained in the
#' \code{estimates} slot.
#' For example r.values for correlations or r.squared for models.
#' @slot df a vector of length 2, containing df1 and df2 corresponding to the
#' F-ratio considered.
#' @slot estimates A matrix containing the estimates for the parameters named in
#' the \code{type} slot.
#' @slot f.values A matrix containing the f-values, for least-squares.
#' @slot t.values A matrix containing the t-values, for correlations.
#' @slot p.values A matrix containing the p.values for the parameters named in
#' the \code{type} slot.
#' @description \code{anansiTale} is the main container that will hold your
#' stats output data coming out of the \code{anansi} pipeline.
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
