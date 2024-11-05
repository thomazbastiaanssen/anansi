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
#' @slot df a vector of length 2, containing df1 and df2 corresponding to the F-ratio considered.
#' @slot F.valules A matrix containing the f-values
#' @slot estimates A matrix containing the estimates for the parameters named in the \code{type} slot.
#' @slot p.values A matrix containing the p.values for the parameters named in the \code{type} slot.
#' @slot q.values A matrix containing the q.values for the parameters named in the \code{type} slot.
#' @description \code{anansiTale} is the main container that will hold your stats output data coming out of the \code{anansi} pipeline.
#'
setClass("anansiTale",
         slots = c(
           subject   = "character",
           type      = "character",
           df        = "numeric",
           estimates = "matrix",
           F.values  = "matrix",
           p.values  = "matrix",
           q.values  = "matrix"
           )
         )



#' An S4 class to contain all \code{anansi} input so that they can easily be extracted.
#'
#' @slot input A list that holds the input data in \code{anansiWeb} format, as well as the \code{groups}, formula and \code{reff} argument(s) if provided.
#' @description \code{anansiInput} is the container that will hold your input data in the \code{anansiYarn} output file coming out of the \code{anansi} pipeline.
#'
setClass("anansiInput",
         slots = c(
           web    = "anansiWeb",
           groups = "vector",
           formula = "formula",
           reff   = "vector"
         )
)

#' An S4 class to contain all \code{anansi} output so that they can easily be extracted.
#'
#' @slot outut A list that holds the output data in \code{anansiWeb} format, as well as the \code{groups} argument if provided.
#' @description \code{anansiOutput} is the container that will hold your output data in the \code{anansiYarn} output file coming out of the \code{anansi} pipeline.
#'
setClass("anansiOutput",
         slots = c(
           cor_results = "list",
           model_results = "list"
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
           input  = "anansiInput",
           output = "anansiOutput"
         )
)

#' Plotting method for \code{anansiYarn} object.
#' @description Makes a lot of histograms, useful to eyeball the p-value distribution.
#' @param x An \code{anansiYarn} object.
#' @importFrom graphics hist par
#' @export
#'
setMethod("plot", "anansiYarn", function(x){
  par(ask=TRUE)
  plotnames  = c(names(x@output@cor_results), names(x@output@model_results))
  plotlist   = list(unlist(x@output@cor_results), unlist(x@output@model_results))
  dictionary = x@input@web@dictionary

  for(p in 1:(length(plotnames))){
    tale = unlist(plotlist)[[plotnames[p]]]

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
         xlim = c(0, 1), xlab = "p.values",
         main = paste("p.values of", tale@subject),
         breaks = 20)

    #plot q.values
    hist(c(tale@q.values[dictionary]),
         xlim = c(0, 1), xlab = "q.values",
         main = paste("q.values of", tale@subject),
         breaks = 20)
  }
  par(ask = FALSE)
  }
  )

#' Show method for \code{anansiYarn} object.
#' @description method to print \code{anansiYarn} object by calling \code{show}.
#' Since anansiYarn objects are typically huge and unwieldy, also gives some tips on how to handle it.
#' @param object An \code{anansiYarn} object.
#' @importFrom utils str
#' @importFrom methods show
#' @export
#'
setMethod("show", "anansiYarn", function(object){
  str(object)
  cat("\nThis is an anansiYarn S4 object. They tend to be very large so we here's a summary instead.
      You can parse it to a data.frame using the spinToWide() & spinToLong() functions,
      or you can manually explore it by using the @ operator. ")
}
)

#' Print method for \code{anansiYarn} object.
#' @description method to print \code{anansiYarn} object by calling \code{print}.
#' Since anansiYarn objects are typically huge and unwieldy, also gives some tips on how to handle it.
#' @param x An \code{anansiYarn} object.
#' @importFrom utils str
#' @export
#'
setMethod("print", "anansiYarn", function(x){
  str(x)
  cat("\nThis is an anansiYarn S4 object. They tend to be very large so we here's a summary instead.
      You can parse it to a data.frame using the spinToWide() & spinToLong() functions,
      or you can manually explore it by using the @ operator. ")
}
)

#' Summary method for \code{anansiYarn} object.
#' @description method to print \code{anansiYarn} object by calling \code{summary}.
#' Since anansiYarn objects are typically huge and unwieldy, also gives some tips on how to handle it.
#' @param object An \code{anansiYarn} object.
#' @importFrom utils str
#' @export
#'
setMethod("summary", "anansiYarn", function(object){
  str(object)
  cat("\nThis is an anansiYarn S4 object. They tend to be very large so we here's a summary instead.
      You can parse it to a data.frame using the spinToWide() & spinToLong() functions,
      or you can manually explore it by using the @ operator. ")
}
)

