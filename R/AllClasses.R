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

#' Plotting method for \code{anansiYarn} object.
#' @description Makes a lot of histograms, useful to eyeball the p-value distribution.
#' @param x An \code{anansiYarn} object.
#' @importFrom graphics hist par
#'
setMethod("plot", "anansiYarn", function(x){
  par(ask=TRUE)
  plotnames  = names(unlist(x@output))
  dictionary = x@input$web@dictionary

  for(p in 1:(length(plotnames))){
    tale = unlist(x@output)[[plotnames[p]]]

    #figure out if we're plotting r.values or r.squared.
    if(tale@type == "r.values"){
    hist(c(tale@estimates[dictionary]),
         xlim = c(-1, 1), xlab = tale@type, main = paste(tale@type, "of", tale@subject))
    }
    if(tale@type == "r.squared"){
      hist(c(tale@estimates[dictionary]),
           xlim = c(0, 1), xlab = tale@type, main = paste(tale@type, "of", tale@subject))
    }
    #plot p.values
    hist(c(tale@p.values[dictionary]),
         xlim = c(0, 1), xlab = tale@type,
         main = paste("p.values of", tale@subject),
         breaks = 20)

    #plot q.values
    hist(c(tale@q.values[dictionary]),
         xlim = c(0, 1), xlab = tale@type,
         main = paste("q.values of", tale@subject),
         breaks = 20)
  }
  par(ask = FALSE)
  }
  )
