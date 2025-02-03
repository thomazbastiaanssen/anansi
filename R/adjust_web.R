#' Partial out covariates using 'lm'.
#' @description Fit a model and return the residuals of the tableY features.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param formula a formula object of type ~ x + y + z where x, y and z are covariates that you are not interested in. Do not include groups or tableX features here.
#' @param metadata a data.frame containing the term(s) mentioned in the formula
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return An anansiWeb object with updated tableY.
#' @export
#'
adjust_web <- function(web, formula, metadata, verbose = T) {
  web@tableY <- pre_partial(web = web, formula = formula, metadata, verbose = verbose)
  return(web)
}

#' Partial out covariates using 'lm'.
#' @description Fit a model and return the residuals of the tableY features. Called by \code{adjust_web()}
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param formula a formula object of type ~ x + y + z where x, y and z are covariates that you are not interested in. Do not include groups or tableX features here.
#' @param metadata a data.frame containing the term(s) mentioned in the formula
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return An anansiWeb object with updated tableY.
#' @importFrom stats lm residuals update.formula
#'
pre_partial <- function(web, formula, metadata, verbose) {
  if (verbose) {
    print(paste("Partialling out the following terms:", paste(formula, collapse = " ")))
  }
  Ypre <- web@tableY
  Ypost <- apply(Ypre, 2, FUN = function(y) {
    temp_meta <- cbind(y_to_meta = y, metadata)
    residuals(lm(formula = update.formula(formula, new = y_to_meta ~ .), data = temp_meta))
  })
  return(Ypost)
}
