#' Extract information from an anansiTale object and parse it into a neat table
#' @param tale An \code{anansiTale} object
#' @param format either \code{wide} or \code{long}. Default is \code{wide}. Toggles the output format.
#' @return A wide format data.frame. The first two columns are the Compounds and Functions.
#' @seealso \code{\link{spinToWide}}\cr \code{\link{spinToLong}}
#'

getAnansiResults <- function(tale, format = "wide"){
  if(!format %in% c("wide", "long")){stop("format should be either `wide` or `long`.")}
  out_df =
  cbind(expand.grid(
  Compounds = row.names(tale@estimates),
  Functions = colnames(tale@estimates), stringsAsFactors = F),
  estimates = c(tale@estimates),
  p.values  = c(tale@p.values),
  q.values  = c(tale@q.values))

  if(format == "wide"){
  names(out_df)[3      ] <- tale@type
  names(out_df)[-c(1:2)] <- paste(tale@subject, names(out_df)[-c(1:2)], sep = "_")
  }
  if(format == "long"){
    out_df$type = paste(tale@subject, sep = "_")
  }

  return(out_df)
}
