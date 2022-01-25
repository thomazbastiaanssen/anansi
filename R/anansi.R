#' Calculate an association network
#' @param web An anansiWeb object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables} for how to weave a web.
#' @param method Correlation method. Method="pearson" is the default value. The alternatives to be passed to cor are "spearman" and "kendall".
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param adjust.method Method to adjust p-values for multiple comparisons. Method="BH" is the default value. See p.adjust in the base R stats package.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return A list of lists containing correlation coefficients, p-values and q-values for all operations.

anansi = function(web, method = "pearson", groups, adjust.method = "BH", verbose = T){

  if(is.null(groups) | is.na(groups)){message("Please be aware that the `groups` argument is missing. \nAnansi will proceed without differential association testing"  )}

  if(verbose){print("Running annotation-based correlations")}
  out_cors  = anansiCorTestByGroup(web = web, method = method, groups = groups, adjust.method = adjust.method, verbose = verbose)

  if(verbose){print("Fitting models for differential correlation testing")}
  out_models = anansiDiffCor(web = web, groups = groups, adjust.method = adjust.method)

  return(list(cor_results = out_cors, model_results = out_models))
}
