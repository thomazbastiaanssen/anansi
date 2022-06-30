#' Run correlations for all interacting metabolites and functions.
#' @description If the \code{groups} argument is suitable, will also run correlation analysis per group. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor()} are "spearman" and "kendall".
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total dataset and per group if applicable.
#' @seealso \code{\link{anansi}}
#'
anansiCorTestByGroup = function(web, method = "pearson", groups, verbose = T){

  #Determine all groups
  all_groups = unique(groups)
  #If there are numbers here, we cannot do cor by group, so we'll substitute a groups called "All" for this part.
  if(!inherits(all_groups, "character")){all_groups = "All" ; groups = rep("All", nrow(web@tableY))}

  #If that's all fine, we can determine the size of the output
  else if(length(all_groups) > 1)
  {all_groups = c("All", all_groups)}

  #Generate container list of suitable length for all results
  out_list = vector(mode = "list", length = length(all_groups))
  names(out_list) = c(all_groups)

  #first run for all groups together
  out_list$All = anansiCorPvalue(web, method = method, groups = rep(T, nrow(web@tableY)), verbose = verbose)
  out_list$All@subject = "All"


  #If verbose, verbalize.
  if(length(all_groups) > 1){

  if(verbose)
  {print(paste("Running correlations for the following groups:",
                    paste(all_groups, collapse = ", ")))}

  #Assess correlations for subsets if applicable.
  #Skip 1 since that's taken by "All".
  for(i in 2:length(all_groups)){
    out_by_group = anansiCorPvalue(web, method = method, groups = groups == all_groups[i], verbose = verbose)
    out_by_group@subject = all_groups[i]

    out_list[[i]] = out_by_group
    }
  }
  #Return results
  return(out_list)
}

#' Compute r-statistics for each featureY-featureX pair in the dictionary.
#' Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor()} are "spearman" and "kendall".
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return An \code{anansiTale} result object.
#' @seealso \code{\link{anansi}} \cr \code{\link{anansiCorTestByGroup}}
#' @importFrom stats pt
#' @importFrom methods new
#'
anansiCorPvalue = function(web, method = "pearson", groups = NULL, verbose = verbose) {
  #Compute correlation coefficients
  r    <- anansiCor(web = web, method = method, groups = groups)

  #Compute t-statistics based on the n and the correlation coefficient
  n    <- web@dictionary
  n[T] <- nrow(web@tableY[groups,])
  t    <- (r*sqrt(n-2))/sqrt(1-r^2)

  #Compute p-values based on t and n.
  p    <- 2*(1 - pt(abs(t),(n-2)))

  #Compute naive adjusted p-values
  q    <- p
  q[web@dictionary] <- NA

  #Collate correlation coefficients, p-values and q-values into an anansiTale
  out  <- new("anansiTale",
              subject    = "All",
              type       = "r.values",
              estimates  = r,
              p.values   = p,
              q.values   = q)
  return(out)
}

#' Compute r-statistics for each featureY-featureX pair in the dictionary.
#' Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor()} are "spearman" and "kendall".
#' @param groups A boolean vector used to select which samples should be included in the correlations.
#' @seealso \code{\link{anansi}} \cr \code{\link{anansiCorTestByGroup}}
#' @return A matrix of r-statistics.
#' @importFrom stats cor
#'
anansiCor = function(web, method = "pearson", groups = NULL){

  #Run correlations on subsections of your data
  merge <- cbind(web@tableY[groups,], web@tableX[groups,])
  cors  <- cor(merge, method = method, use = "pairwise.complete.obs")
  cors_bipartite <- cors[1:ncol(web@tableY), (ncol(web@tableY)+1):ncol(cors)]
  #set non-canonical correlations to zero using the binary adjacency matrix.
  return(cors_bipartite * web@dictionary)
}
