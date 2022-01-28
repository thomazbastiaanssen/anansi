#' Run correlations for all interacting metabolites and functions.
#' @description If the \code{groups} argument is suitable, will also run correlation analysis per group. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor()} are "spearman" and "kendall".
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param adjust.method Method to adjust p-values for multiple comparisons. \code{adjust.method = "BH"} is the default value. See \code{p.adjust()} in the base R \code{stats} package.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total dataset and per group if applicable.
#' @seealso \code{\link{anansi}}
#'
anansiCorTestByGroup = function(web, method = "pearson", groups, adjust.method = adjust.method, verbose = T){
  #Determine all groups
  all_groups = unique(groups)

  #Assess correlations for entire data set
  all_out = anansiCorPvalue(web, method = method, adjust.method = adjust.method)
  all_out@subject = "All"
  out_list = list(All = all_out)

  #If groups argument is suitable, run subset analysis
  if(inherits(groups, "character") & all(table(groups) > 3) & length(all_groups > 1)){

  #If verbose, verbalize.
  if(verbose){print(paste("Running correlations for the following groups:",
                          paste(all_groups, collapse = ", "),
                          "and all together"))}

  #Assess correlations for subsets if applicable.
  for(i in 1:length(all_groups)){
    out_by_group = anansiCorPvalue(web, method = method, groups = groups == all_groups[i], adjust.method = adjust.method)
    out_by_group@subject = all_groups[i]

    out_list[[i+1]] = out_by_group
  }
  names(out_list)[-1] <- all_groups
  }
  #Return results
  return(out_list)
}

#' Compute r-statistics for each compound-function pair in the dictionary.
#' Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor()} are "spearman" and "kendall".
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param adjust.method Method to adjust p-values for multiple comparisons. \code{adjust.method = "BH"} is the default value. See \code{p.adjust()} in the base R \code{stats} package.
#' @return An \code{anansiTale} result object.
#' @seealso \code{\link{anansi}} \cr \code{\link{anansiCorTestByGroup}}
#' @importFrom stats pt p.adjust
#' @importFrom methods new
#'
anansiCorPvalue = function(web, method = "pearson", groups = NULL, adjust.method = adjust.method) {
  #Compute correlation coefficients
  r    <- anansiCor(web = web, method = method, groups = groups)

  if(is.null(groups) | any(is.na(groups))){groups = TRUE}
  #Compute t-statistics based on the n and the correlation coefficient
  n    <- web@dictionary
  n[T] <- nrow(web@tableY[groups,])
  t    <- (r*sqrt(n-2))/sqrt(1-r^2)

  #Compute p-values based on t and n.
  p    <- 2*(1 - pt(abs(t),(n-2)))

  #Compute naive adjusted p-values
  q    <- p
  q[web@dictionary] <- p.adjust(p[web@dictionary], method = adjust.method)

  #Collate correlation coefficients, p-values and q-values into an anansiTale
  out  <- new("anansiTale",
              subject    = "All",
              type       = "r.values",
              estimates  = r,
              p.values   = p,
              q.values   = q)
  return(out)
}

#' Compute r-statistics for each compound-function pair in the dictionary.
#' Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor()} are "spearman" and "kendall".
#' @param groups A boolean vector used to select which samples should be included in the correlaitons.
#' @seealso \code{\link{anansi}} \cr \code{\link{anansiCorTestByGroup}}
#' @return A matrix of r-statistics.
#' @importFrom stats cor
#'
anansiCor = function(web, method = "pearson", groups = NULL){
  #Run correlations on subsections of your data
  if(is.null(groups) | any(is.na(groups))){groups = TRUE}
  merge <- cbind(web@tableY[groups,], web@tableX[groups,])
  cors  <- cor(merge, method = method)
  cors_bipartite <- cors[1:ncol(web@tableY), (ncol(web@tableY)+1):ncol(cors)]
  #set non-canonical correlations to zero using the binary adjacency matrix.
  return(cors_bipartite * web@dictionary)
}
