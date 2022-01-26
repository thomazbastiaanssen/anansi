#' Calculate an association network
#' @param web An anansiWeb object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables} for how to weave a web.
#' @param method Correlation method. Method="pearson" is the default value. The alternatives to be passed to cor are "spearman" and "kendall".
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param adjust.method Method to adjust p-values for multiple comparisons. Method="BH" is the default value. See p.adjust in the base R stats package.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @param diff_cor A boolean. Toggles whether to compute differential correlations. Default is True.
#' @return A list of lists containing correlation coefficients, p-values and q-values for all operations.

anansi = function(web, method = "pearson", groups, adjust.method = "BH", verbose = T, diff_cor = T){
  if(inherits(groups, "character")){
    if(!all(table(groups) > 3)) {
      warning("The `groups` argument is categorical, but not all groups have at least three observations.
              Correlations per group cannot be done. Please check your groups. ")
      diff_cor = F
      groups = "All"
    }
  }

  if(diff_cor & (is.null(groups) | any(is.na(groups)))){
    message("Please be aware that the `groups` argument is missing or contains NAs. \n
            Anansi will proceed without differential association testing")
    diff_cor = F}

  if(verbose){print("Running annotation-based correlations")}
  outlist = list(cor_results = anansiCorTestByGroup(web = web,
                                                    method = method,
                                                    groups = groups,
                                                    adjust.method = adjust.method,
                                                    verbose = verbose))

  if(diff_cor){
  if(verbose){print("Fitting models for differential correlation testing")}
    outlist[["model_results"]] = anansiDiffCor(web = web, groups = groups, adjust.method = adjust.method)
  }
  return(outlist)
}
