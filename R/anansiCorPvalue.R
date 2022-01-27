#' Compute r-statistics for each compound-function pair in the dictionary.
#' Typically, the main \code{anansi} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables} for how to weave a web.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor} are "spearman" and "kendall".
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param adjust.method Method to adjust p-values for multiple comparisons. \code{adjust.method = "BH"} is the default value. See \code{p.adjust} in the base R \code{stats} package.
#' @return An \code{anansiTale} result object.
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
             type       = "correlation",
             estimates  = r,
             p.values   = p,
             q.values   = q)
  return(out)
}
