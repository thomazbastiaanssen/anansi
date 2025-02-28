#' Handle FDR methods for anansi.
#' @param results The results table containing p-values to be adjusted.
#' @param method The p-value adjustment method. See ?p.adjust.
#' @return The `results` table, extended with adjusted p.values.
#' @importFrom stats p.adjust
#'
anansi.p.adjust <- function(results, method){
  p.cols <- grep("p.values", colnames(results))
  q.cols <- apply(results[,p.cols, drop = FALSE], 2, p.adjust, method)
  colnames(q.cols) <- gsub(
    "_p.values", "_q.values", x = colnames(q.cols), fixed = TRUE)
  cbind.data.frame(results, q.cols)
}
