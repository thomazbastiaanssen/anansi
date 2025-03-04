#' Extract information from an anansiTale object and parse it into a neat table
#' @param tale An \code{anansiTale} object
#' @param dic A dictionary.
#' @noRd
#' @return A wide format data.frame with summary statistics by feature pair.
#'
frame.tale <- function(tale, dic){
  switch(tale@type,
         "r.values"  = frame.tale.cor(tale, dic),
         "r.squared" = frame.tale.ols(tale, dic))
}

#' @noRd
#'
frame.tale.cor <- function(tale, dic){
  out.df <- data.frame(
    r.values = tale@estimates[dic],
    t.values = tale@t.values[dic],
    p.values = tale@p.values[dic]
    #q.values =  tale@q.values[dic]
    )
  colnames(out.df) <- paste(tale@subject, colnames(out.df), sep = "_")
  return(out.df)
}

#' @noRd
#'
frame.tale.ols <- function(tale, dic){
  out.df <- data.frame(
    r.squared = tale@estimates[dic],
    f.values =  tale@f.values[dic],
    p.values =  tale@p.values[dic]
    #q.values =  tale@q.values[dic]
    )
  colnames(out.df) <- paste(tale@subject, colnames(out.df), sep = "_")
  return(out.df)
}

#' Shape list of `anansiTale` results into a data.frame.
#' @noRd
#'
result.df <- function(out.list, dic) {
  feature_labs <- expand.grid(
         feature_Y = row.names(dic),
         feature_X = colnames(dic), stringsAsFactors = FALSE
     )[dic,]

  df.list <- c(feature_labs, lapply(out.list, frame.tale, dic))
  do.call(what = "cbind.data.frame", args = df.list, quote = TRUE)
}

#' Handle FDR methods for anansi.
#' @param results The results table containing p-values to be adjusted.
#' @param method The p-value adjustment method. See ?p.adjust.
#' @return The `results` table, extended with adjusted p.values.
#' @importFrom stats p.adjust
#' @noRd
#'
anansi.p.adjust <- function(results, method){
    p.cols <- grep("p.values", colnames(results))
    q.cols <- apply(results[,p.cols, drop = FALSE], 2, p.adjust, method)
    colnames(q.cols) <- gsub(
        "_p.values", "_q.values", x = colnames(q.cols), fixed = TRUE)
    cbind.data.frame(results, q.cols)
}