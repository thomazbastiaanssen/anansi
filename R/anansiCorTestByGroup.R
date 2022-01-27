#' Run correlations for all interacting metabolites and functions.
#' If the \code{groups} argument is suitable, will also run correlation analysis per group. Typically, the main \code{anansi} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables} for how to weave a web.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor} are "spearman" and "kendall".
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param adjust.method Method to adjust p-values for multiple comparisons. \code{adjust.method = "BH"} is the default value. See \code{p.adjust} in the base R \code{stats} package.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total dataset and per group if applicable.
#'
anansiCorTestByGroup = function(web, method = "pearson", groups, adjust.method = adjust.method, verbose = T){
  #Determine all groups
  all_groups = unique(groups)

  #Assess correlations for entire data set
  all_out = anansiCorPvalue(web, method = method, adjust.method = adjust.method)
  all_out@type = paste(all_out@type, "All", sep = "_")
  out_list = list(All = all_out)

  #If groups argument is suitable, run subset analysis


  ######################
  ###################### fix me!!
  ######################
  if(length(unique(groups)) < floor(length(groups)/3)){

  #If verbose, verbalize.
  if(verbose){print(paste("Running correlations for the following groups:",
                          paste(all_groups, collapse = ", "),
                          "and all together"))}

  #Assess correlations for subsets if applicable.
  for(i in 1:length(all_groups)){
    out_by_group = anansiCorPvalue(web, method = method, groups = groups == all_groups[i], adjust.method = adjust.method)
    out_by_group@type = paste(out_by_group@type, all_groups[i], sep = "_")

    out_list[[i+1]] = out_by_group
  }
  names(out_list)[-1] <- all_groups
  }
  #Return results
  return(out_list)
}
