#' Take anansi output and wrangle it all to a wide data.frame format.
#' @param anansi_output The output of the main anansi function.
#' @param prune Boolean, default is TRUE. Toggles whether to take out the non-canonical associations.
#' @param translate Boolean, default is TRUE. Toggles whether to translate the names of Compounds and Functions to human readable names.
#' @param Y_translation data.frame, a lookup table with compound names as the first column and human readable names as the second. See \code{cpd_translation} for an example file.
#' @param X_translation data.frame, a lookup table with function names as the first column and human readable names as the second. See \code{KO_translation} for an example file.
#' @return a wide format data.frame
#' @export
#'
spinToWide <- function(anansi_output, prune = T, translate = T, Y_translation = anansi::cpd_translation, X_translation = anansi::KO_translation){
  #First flatten all types  of results (cor, model, etc) and create a list of individual wide data.frames.
  flat_list = lapply(unlist(list(anansi_output@output@cor_results,
                                 anansi_output@output@model_results)), getAnansiResults)


  #merge all data.frames in the list while keeping the row order.
  wide_df   = Reduce(function(x, y) merge(x, y, sort = F), flat_list)

  if(translate){
    wide_df = anansiTranslate(wide_df, Y_translation = Y_translation, X_translation = X_translation)
  }

  if(prune){
    #If true, remove all non-canonical interactions.
    wide_df = wide_df[c(anansi_output@input@web@dictionary),]
  }

  return(wide_df)
}


#' Take anansi output and wrangle it all to a long data.frame format.
#' @param anansi_output an \code{anansiYarn} object. The output of the main anansi function.
#' @param prune Boolean, default is TRUE. Toggles whether to take out the non-canonical associations.
#' @param translate Boolean, default is TRUE. Toggles whether to translate the names of Compounds and Functions to human readable names.
#' @param Y_translation data.frame, a lookup table with compound names as the first column and human readable names as the second. See \code{cpd_translation} for an example file.
#' @param X_translation data.frame, a lookup table with function names as the first column and human readable names as the second. See \code{KO_translation} for an example file.
#' @return a long format data.frame intended to be compatible with \code{ggplot2}
#' @export
#'
spinToLong <- function(anansi_output, prune = T, translate = T, Y_translation = anansi::cpd_translation, X_translation = anansi::KO_translation){

  #Figure out how many groups were present in the correlations.
  n_cors = length(anansi_output@output@cor_results)

  #Flatten all types  of the correlations and create a list of individual wide data.frames.
  flat_cor_list = lapply(anansi_output@output@cor_results, getAnansiResults, format = "long")

  #Merge all individual correlation results into a single long data.frame
  long_out  = Reduce(rbind, flat_cor_list)
  colnames(long_out)[3] = "r.values"

  if(length(anansi_output@output@model_results) > 0){
    #Make a flat list of the model results in wide format
    flat_model_list = lapply(anansi_output@output@model_results, getAnansiResults, format = "wide")

    #Merge all model results in a single wide data.frame
    wide_model_df   = Reduce(function(x, y) merge(x, y, sort = F), flat_model_list)

    #Stack the model results so that they have the same number of rows as the correlation results
    wide_model_df   = do.call("rbind", replicate(n_cors, wide_model_df, simplify = FALSE))

    #Combine the results into a single data.frame
    long_out = cbind(long_out, wide_model_df[,-c(1, 2)])
  }

  if(translate){
    long_out = anansiTranslate(long_out, Y_translation = Y_translation, X_translation = X_translation)
  }

  if(prune){
    #If true, remove all non-canonical interactions.
    long_out = long_out[c(anansi_output@input@web@dictionary),]

  }

  return(long_out)
}

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


#' Translate Compounds and Functions column to human readable names
#' @description helper funciton for the spinToLong and \code{spinToWide} functions.
#' @param x a table with a Compounds and Functions column.
#' @param Y_translation data.frame, a lookup table with compound names as the first column and human readable names as the second. See \code{cpd_translation} for an example file.
#' @param X_translation data.frame, a lookup table with function names as the first column and human readable names as the second. See \code{KO_translation} for an example file.
#' @return an expanded table with a Compounds and Functions column.
#'
anansiTranslate <- function(x, Y_translation = Y_translation, X_translation = X_translation){

  relevant_cpd = Y_translation[Y_translation[,1] %in% x$Compounds,]
  x$Compounds = rep(paste(relevant_cpd[,1],
                          gsub(";.*", "", relevant_cpd[,2]), sep = ": "),
                    times = length(unique(x$Functions)))

  relevant_fun = X_translation[X_translation[,1] %in% x$Functions,]
  x$Functions = rep(paste(relevant_fun[,1],
                          relevant_fun[,2], sep = ": "),
                    each = length(unique(x$Compounds)))
  return(x)
}
