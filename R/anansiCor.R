#' Compute r-statistics for each compound-function pair in the dictionary.
#' Typically, the main \code{anansi} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables} for how to weave a web.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor} are "spearman" and "kendall".
#' @return A matrix of r-statistics.
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
