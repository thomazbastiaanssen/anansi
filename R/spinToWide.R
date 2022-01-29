#' Take anansi output and wrangle it all to a wide data.frame format.
#' @param anansi_output The output of the main anansi function.
#' @param prune Boolean, default is TRUE. Toggles whether to take out the non-canonical associations.
#' @return a wide format data.frame
#' @export

spinToWide <- function(anansi_output, prune = T){
  #First flatten all types  of results (cor, model, etc) and create a list of individual wide data.frames.
  flat_list = lapply(unlist(anansi_output@output), getAnansiResults)


  #merge all data.frames in the list while keeping the row order.
  wide_df   = Reduce(function(x, y) merge(x, y, sort = F), flat_list)

  if(prune){
    #If true, remove all non-canonical interactions.
    wide_df = wide_df[c(anansi_output@input$web@dictionary),]
  }

  return(wide_df)
}
