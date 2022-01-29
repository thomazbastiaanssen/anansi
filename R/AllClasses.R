#' An S4 class to contain all metabolomics and functional input data as well as a dictionary to link them.
#'
#' @slot tableY A matrix of metabolomics data. Rows are samples and columns are features.
#' @slot tableX A matrix of functional data. Rows are samples and columns are features.
#' @slot dictionary A binary adjacency matrix. Typically generated using the \code{weaveWebFromTables()} function.
#' @description anansiWeb is the main container that will hold your input data thoughout the \code{anansi} pipeline.
#'
setClass("anansiWeb",
         slots = c(
           tableY     = "matrix",
           tableX     = "matrix",
           dictionary = "matrix"
           )
         )

#' An S4 class to contain all \code{anansi} stats results so that they can easily be extracted.
#'
#' @slot subject A character that describes the data that was queried.
#' @slot type A character that describes type of parameter contained in the \code{estimates} slot.
#' For example r.values for correlations or r.squared for models.
#' @slot estimates A data.frame containing the estimates for the parameters named in the \code{type} slot.
#' @slot p.values A data.frame containing the p.values for the parameters named in the \code{type} slot.
#' @slot q.values A data.frame containing the q.values for the parameters named in the \code{type} slot.
#' @description \code{anansiTale} is the main container that will hold your stats output data coming out of the \code{anansi} pipeline.
#'
setClass("anansiTale",
         slots = c(
           subject   = "character",
           type      = "character",
           estimates = "matrix",
           p.values  = "matrix",
           q.values  = "matrix"
           )
         )

#' An S4 class to contain all \code{anansi} input and output so that they can easily be extracted.
#'
#' @slot input A list that holds the input data in \code{anansiWeb} format, as well as the \code{groups} argument if provided.
#' @slot output A list that holds the output in \code{anansiTale} format.
#' @description \code{anansiYarn} is the main container that will hold your output data coming out of the \code{anansi} pipeline.
#' The \code{spinToWide()} and \code{spinToLong()} functions can be used to extract result tables in wide and long format, respectively.
#'
setClass("anansiYarn",
         slots = c(
           input  = "list",
           output = "list"
         )
)
