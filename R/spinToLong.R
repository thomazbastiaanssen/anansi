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
  #Figure out how many types of results were computed (cor, model, etc). be sure to exclude the dictionary here
  res_types = length(anansi_output@output)
  res_names = names(anansi_output@output)[1:res_types]

  #Also, figure out how many groups were present in the correlations.
  n_cors = length(anansi_output@output$cor_results)

  #Flatten all types  of the correlations and create a list of individual wide data.frames.
  flat_cor_list = lapply(anansi_output@output$cor_results, getAnansiResults, format = "long")

  #Merge all individual correlation results into a single long data.frame
  long_out  = Reduce(rbind, flat_cor_list)
  colnames(long_out)[3] = "r.values"

  if("model_results" %in% res_names){
    #Make a flat list of the model results in wide format
    flat_model_list = lapply(anansi_output@output$model_results, getAnansiResults, format = "wide")

    #Merge all model results in a single wide data.frame
    wide_model_df   = Reduce(function(x, y) merge(x, y, sort = F), flat_model_list)

    #Stack the model results so that they have the same number of rows as the correlation results
    wide_model_df   = do.call("rbind", replicate(n_cors, wide_model_df, simplify = FALSE))

    #Combine the results into a single data.frame
    long_out = cbind(long_out, wide_model_df[,-c(1, 2)])
  }

  if(translate){
    relevant_cpd = Y_translation[Y_translation[,1] %in% long_out$Compounds,]
    long_out$Compounds = rep(paste(relevant_cpd[,1],
                                  gsub(";.*", "", relevant_cpd[,2]), sep = ": "),
                            times = length(unique(long_out$Functions)))

    relevant_fun = X_translation[X_translation[,1] %in% long_out$Functions,]
    long_out$Functions = rep(paste(relevant_fun[,1],
                                  relevant_fun[,2], sep = ": "),
                            each = length(unique(long_out$Compounds)))
  }

  if(prune){
    #If true, remove all non-canonical interactions.
    long_out = long_out[c(anansi_output@input$web@dictionary),]

  }

  return(long_out)
}
