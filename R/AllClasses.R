#' An S4 class to contain all metabolomics and functional input data as well as
#' a dictionary to link them.
#' @slot tableY A matrix of metabolomics data. Rows are samples and columns are
#' features.
#' @slot tableX A Matrix of functional data. Rows are samples and columns are
#' features.
#' @slot dictionary A binary adjacency matrix. Typically generated using the
#' \code{weaveWebFromTables()} function.
#' @description \code{anansiWeb} is the main container that will hold your input
#' data throughout the \code{anansi} pipeline.
#'
setClass("anansiWeb",
  slots = c(
    tableY     = "ANY",
    tableX     = "ANY",
    dictionary = "ANY"
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
#' @slot f.valules A matrix containing the f-values, for least-squares.
#' @slot t.valules A matrix containing the t-values, for correlations.
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
    estimates = "ANY",
    f.values  = "ANY",
    t.values  = "ANY",
    p.values  = "ANY"
  )
)

#' An S4 class to contain all \code{anansi} input so that they can easily be
#' extracted.
#'
#' @slot input A list that holds the input data in \code{anansiWeb} format, as
#' well as the \code{int.terms} and corresponding \code{groups}, the equivalent
#'\code{lm.formula} and \code{error.term} argument(s) if provided.
#' @description \code{anansiInput} is the container that will hold your input
#' data in the \code{anansiYarn} output file coming out of the \code{anansi}
#' pipeline.
#'
setClass("anansiInput",
  slots = c(
    web = "anansiWeb",
    lm.formula = "formula",
    error.term = "ANY",
    int.terms = "vector",
    groups = "ANY",
    n.grps = "vector",
    group.id = "vector",
    metadata = "ANY"
  )
)
