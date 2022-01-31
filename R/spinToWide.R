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
  flat_list = lapply(unlist(anansi_output@output), getAnansiResults)


  #merge all data.frames in the list while keeping the row order.
  wide_df   = Reduce(function(x, y) merge(x, y, sort = F), flat_list)

  if(translate){
    relevant_cpd = Y_translation[Y_translation[,1] %in% wide_df$Compounds,]
    wide_df$Compounds = rep(paste(relevant_cpd[,1],
                            gsub(";.*", "", relevant_cpd[,2]), sep = ": "),
                            times = length(unique(wide_df$Functions)))

    relevant_fun = X_translation[X_translation[,1] %in% wide_df$Functions,]
    wide_df$Functions = rep(paste(relevant_fun[,1],
                                  relevant_fun[,2], sep = ": "),
                            each = length(unique(wide_df$Compounds)))
  }

  if(prune){
    #If true, remove all non-canonical interactions.
    wide_df = wide_df[c(anansi_output@input$web@dictionary),]
  }

  return(wide_df)
}
