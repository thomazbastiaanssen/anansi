#'Wrangle anansi dictionary list into binary adjacency matrix
#' @param tableY A matrix containing metabolites of interest. Rows should be samples and columns should be features.
#' @param dictionary A list that has compound names as entries. For general use, we recommend using the one provided in this package.
#' @return a binary adjacency matrix with compounds form tableY as rows and functions from tableX as columns.
#'
makeAdjacencyMatrixFromList <- function(tableY, dict_list){
  #Prune list to only contain metabolites in tableY
  dict_pruned = dict_list[names(dict_list) %in% colnames(tableY)]
  compnames   = names(dict_pruned)
  funcnames   = sort(unique(unlist(dict_pruned)))

  #create an empty adjacency matrix
  dict_out = matrix(nrow = length(compnames),
                    ncol = length(funcnames),
                    dimnames = list(compnames, funcnames),
                    data = F)
  #Fill in the canonical associations by row
  for(comp in 1:length(dict_pruned)){
    dict_out[comp,dict_pruned[[comp]]] <- T
  }
  #Return the resulting matrix
  return(dict_out)
}
