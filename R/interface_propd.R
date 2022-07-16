#' Run propr for all interacting features such as metabolites or functions.
#' Wraps around the propr function.
#' @description If the \code{groups} argument is suitable, will also run proportionality analysis per group. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param groups A categorical or continuous value necessary for differential proportionality Typically a state or treatment score.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total dataset and per group if applicable.
#' @importFrom propr propr
#' @seealso \code{\link{anansi}}
#'
wrap_propr = function(web, groups, verbose = T){
  warning("propr not yet supported.")
  break
  pr <- propr(web@tableY,     # rows as samples, like it should be
              metric = "rho", # or "phi", "phs", "cor", "vlr"
              ivar = NA, )    # used by updateCutoffs
}

#' Run differential proportionality analysis using propr for all interacting features such as metabolites or functions.
#' Wraps around the propd function.
#' @description Can either take continuous or categorical data for \code{groups}. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param groups A categorical or continuous value necessary for differential proportionality Typically a state or treatment score.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent proportionality and one for disjointed proportionality.
#' @importFrom propr propd
#' @importFrom future.apply future_apply
#'
wrap_propd = function(web, groups, verbose = T){
  warning("propr not yet supported.")
  break
  #perhaps softmax to undo CLR? maybe for both?
  pd <- propd(counts = web@tableY,
              group  = groups, # a vector of 2 or more groups
              alpha = NA, # whether to handle zeros
              weighted = TRUE, # whether to weigh log-ratios
              p = 100) # used by updateCutoffs
}
