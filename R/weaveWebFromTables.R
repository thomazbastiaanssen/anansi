#' Create an anansiWeb object from two 'omics tables and a dictionary
#' @description This fucntion will take two dables of 'omics data, typically metabolomics and functional microbiome data. It will also take a dictionary list, which is provided in this package.
#' @param tableY A table containing metabolites of interest. Rows should be samples and columns should be features.
#' @param tableX A table containing functions of interest. Rows should be samples and columns should be features.
#' @param dictionary A list that has compound names from tableY as names. Default is the dictionary provided in this package.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.

#' For general use, we recommend sticking to that one. You can access the dictionary like this: \code{data(dictionary)}
#' @return an \code{anansiWeb} object. Web is used as input for most of the main workflow of anansi.
#' @export
weaveWebFromTables = function(tableY, tableX, dictionary = anansi::anansi_dic, verbose = T){
  #For conventional use, table Y should be metabolites and table X functions.

  # The two 'table' matrices MUST have row
  # and column names that are unique, and
  # look like the following:
  #
  #               f1  f2  f3  f4  f5  f6
  #   sample1      0   0   2   0   0   1
  #   sample2     20   8  12   5  19  26
  #   sample3      3   0   2   0   0   0
  #       ... many more rows ...
  #


  #create binary adjacency matrix first
  dictionary = makeAdjacencyMatrixFromList(tableY = tableY, dict_list = dictionary)

  #Check if input tables have the same names and the same length.
  if(!identical(row.names(tableY), row.names(tableX)))
    {warning("The row names of tableY and tableX do not correspond. Please make sure they are in the same order.")}
  if(nrow(tableY) != nrow(tableX))
    {stop("tableY and tableX do not have the same number of rows/observations.")}

  available_tableY = sort(intersect(colnames(tableY), rownames(dictionary)))
  available_tableX = sort(intersect(colnames(tableX), colnames(dictionary)))

  if(verbose){
    print(paste(length(available_tableY), "were matched between table 1 and the columns of the adjacency matrix"))
    print(paste(length(available_tableX), "were matched between table 2 and the rows of the adjacency matrix"))
  }
  #Select the relevant part of the adjacency matrix
  dictionary = dictionary[available_tableY, available_tableX]

  #Ensure there are no metabolites or functions that never interact
  dictionary = dictionary[rowSums(dictionary) > 0,colSums(dictionary) > 0]

  #use the dictionary to clean input tables
  tableY = tableY[,rownames(dictionary)]
  tableX = tableX[,colnames(dictionary)]

  #Return an anansiWeb object with three slots: typically metabolites, functions and adjacency matrix
  return(new("anansiWeb",
             tableY     = as.matrix(tableY),
             tableX     = as.matrix(tableX),
             dictionary = as.matrix(dictionary)))
}

#'Wrangle anansi dictionary list into binary adjacency matrix
#' @description Takes the anansi dictionary in list format and wrangles it into a biary adjacency matrix based on which compounds are present in both the dictionary and  \code{tableY}.
#' For general use, should probably not be called directly, but rather through \code{\link{weaveWebFromTables}}.
#' @seealso \code{\link{weaveWebFromTables}}
#' @param tableY A matrix containing metabolites of interest. Rows should be samples and columns should be features.
#' @param dict_list A list that has compound names as entries. For general use, we recommend using the one provided in this package.
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
