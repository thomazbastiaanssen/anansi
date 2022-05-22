#' Create an anansiWeb object from two 'omics tables and a dictionary
#' @description This function will take two tables of 'omics data, for example metabolomics and functional microbiome data. It will also take a dictionary list, which is provided in this package.
#' @param tableY A table containing features of interest. Rows should be samples and columns should be features. The Y and X refer to the position of the features in a formula: Y ~ X.
#' @param tableX A table containing features of interest. Rows should be samples and columns should be features. The Y and X refer to the position of the features in a formula: Y ~ X.
#' @param dictionary A list that has feature names from tableY as names. Default is the dictionary provided in this package.
#' @param mode A character vector. Can be "interaction" or "membership". Toggles whether to link two datasets based on their interactions or based on shared group membership.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @param prune A boolean. Toggles whether to prune particularly large groups.
#' @param max_sds A numeric. Only relevant for prune == TRUE. How many SDs larger than the median a group can be before it is pruned.
#' For general use, we recommend sticking to that one. You can access the dictionary like this: \code{data(dictionary)}
#' @return an \code{anansiWeb} object. Web is used as input for most of the main workflow of anansi.
#' @export
weaveWebFromTables = function(tableY, tableX, dictionary = anansi::anansi_dic, verbose = T, mode = "interaction", prune = F, max_sds = 3){

  stopifnot("the mode argument needs to be interaction or membership." = mode %in% c("interaction", "membership"))
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

  if(prune){
    dictionary = prune_dictionary_for_exclusivity(dict_list = dictionary,
                                                  max_sds = max_sds, verbose = verbose)
  }

  #create binary adjacency matrix first
  dictionary = makeAdjacencyMatrix(tableY = tableY, dict_list = dictionary,
                                   verbose = verbose, mode = mode)

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

  #Ensure there are no features that never interact
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

#' Starting point to wrangle anansi dictionary list into binary adjacency matrix
#' @description Takes the anansi dictionary in list format and wrangles it into a binary adjacency matrix.
#' @seealso \code{\link{weaveWebFromTables}}
#' @param tableY A matrix containing y features of interest. Rows should be samples and columns should be features. Only used for mode == "interaction".
#' @param dict_list A list that has tableY names as entries in the case of mode == "interaction" and group names in the case of mode == "membership". For general use, we recommend using the one provided in this package.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @param mode A character vector. Can be "interaction" or "membership". Toggles whether to link two datasets based on their interactions or based on shared group membership.
#' @return a binary adjacency matrix with features from tableY as rows and features from tableX as columns.
#'
makeAdjacencyMatrix <- function(dict_list, tableY = NULL, verbose = T, mode = mode){
  if(mode == "interaction"){
    if(verbose){print("Operating in interaction mode")}
    makeAdjacencyMatrixFromList(tableY = tableY, dict_list = dict_list, mode = mode)
  }
  else if(mode == "membership"){
    if(verbose){print("Operating in membership mode")}
    makeAdjacencyMatrixFromGroupMemberList(dict_list = dict_list, mode = mode)
  }
}

#' Wrangle anansi dictionary list into binary adjacency matrix
#' @description Takes the anansi dictionary in list format and wrangles it into a binary adjacency matrix based on which tableY features are present in both the dictionary and  \code{tableY}.
#' For general use, should probably not be called directly, but rather through \code{\link{weaveWebFromTables}}.
#' @seealso \code{\link{weaveWebFromTables}}
#' @param tableY A matrix containing features of interest. Rows should be samples and columns should be features.
#' @param dict_list A list that has tableY names as entries in the case of mode == "interaction" and groupnames in the case of mode == "membership". For general use, we recommend using the one provided in this package.
#' @param mode A character vector. Can be "interaction" or "membership". Toggles whether to link two datasets based on their interactions or based on shared group membership.
#' @return a binary adjacency matrix with the group members on both the rows and columns.
#'
makeAdjacencyMatrixFromList <- function(tableY = NULL, dict_list, mode = "interaction"){

  #If in interaction mode, prune list to only contain features in tableY
  if(mode == "interaction"){
    dict_list = dict_list[names(dict_list) %in% colnames(tableY)]
  }
  ynames = names(dict_list)
  xnames = sort(unique(unlist(dict_list)))

  #create an empty adjacency matrix
  dict_out = matrix(nrow = length(ynames),
                    ncol = length(xnames),
                    dimnames = list(ynames, xnames), data = F)
  #Fill in the canonical associations by row
  for (y in 1:length(dict_list)) {
    dict_out[y, dict_list[[y]]] <- T
  }
  #Return the resulting matrix
  return(dict_out)
}


#' Wrangle anansi dictionary list into binary adjacency matrix
#' @description Takes the anansi dictionary in list format and wrangles it into a binary adjacency matrix based on group membership.
#' For general use, should probably not be called directly, but rather through \code{\link{weaveWebFromTables}}.
#' @seealso \code{\link{weaveWebFromTables}}
#' @param dict_list A list that has tableY names as entries in the case of mode == "interaction" and groupnames in the case of mode == "membership". For general use, we recommend using the one provided in this package.
#' @param mode A character vector. Can be "interaction" or "membership". Toggles whether to link two datasets based on their interactions or based on shared group membership.
#' @return a binary adjacency matrix with the group members on both the rows and columns.
#'
makeAdjacencyMatrixFromGroupMemberList <- function(dict_list, mode = mode){

  basic_list = makeAdjacencyMatrixFromList(tableY = NULL, dict_list = dict_list, mode = mode)
  #Populate the empty adjacency matrix with TRUE between features that share a membership
  #Thanks to Benjamin Valderrama for the suggestion to use crosspod here!
  dict_out = crossprod(basic_list) >= 1

  return(dict_out)

}

#' Kick out particularly large groups before wrangling the dictionary list into a binary adjacency matrix.
#' @description Takes the anansi dictionary in list format and prunes the larges groups based on sd from the median.
#' @param dict_list A list that has tableY names as entries in the case of mode == "interaction" and group names in the case of mode == "membership". For general use, we recommend using the one provided in this package.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @param max_sds a numeric, how many standard deviantions larger than the group median are groups allowed to be?
#' @return a pruned anansi dictionary list object.
#' @importFrom stats sd
#'
prune_dictionary_for_exclusivity <- function(dict_list, max_sds = 3, verbose = T){

  mem_sizes = unlist(lapply(dict_list, length))
  discard   = mem_sizes > (median(mem_sizes) + (max_sds * sd(mem_sizes)))

  if(verbose)
  {print(paste(sum(discard), "groups were", max_sds, "sds larger than the median group size and were kicked out. "))
    print(paste("These include", paste(names(dict_list[discard]), collapse = ", ")))}

  return(dict_list[!discard])

}
