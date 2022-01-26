#' Extract information from an anansiTale object and parse it into a neat table
#' @param tale An anansiTale object
#' @return A wide format data.frame. The first two columns are the Compounds and Fucntions.
#'

getAnansiResults <- function(tale){

  out_df =
  cbind(expand.grid(
  Compounds = row.names(tale@estimates),
  Functions = colnames(tale@estimates), stringsAsFactors = F),
  estimates = c(tale@estimates),
  p.values  = c(tale@p.values),
  q.values  = c(tale@q.values))

  names(out_df)[-c(1:2)] <- paste(tale@type, names(out_df)[-c(1:2)], sep = "_")

  return(out_df)
}
