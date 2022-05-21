#' Create an anansiWeb object from two 'omics tables and a dictionary
#' @description This fucntion will take two dables of 'omics data, typically metabolomics and functional microbiome data. It will also take a dictionary list, which is provided in this package.
#' @param tableY A table containing metabolites of interest. Rows should be samples and columns should be features.
#' @param tableX A table containing functions of interest. Rows should be samples and columns should be features.
#' @param dictionary A list that has compound names from tableY as names. Default is the dictionary provided in this package.
#' @param max_factor A proportion (0-1). The maximum proportion of the total dataset that can be a member of a group. Only used for membership adjacency matrices.
#' @param clean A boolean, whether to clear temporary files during the function. Only used for membership adjacency matrices
#' @param mode A character vector. Can be "interaction" or "membership". Toggles whether to link two datasets based on their interactions or based on shared group membership.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.

#' For general use, we recommend sticking to that one. You can access the dictionary like this: \code{data(dictionary)}
#' @return an \code{anansiWeb} object. Web is used as input for most of the main workflow of anansi.
#' @export
weaveWebFromTables = function(tableY, tableX, dictionary = anansi::anansi_dic, max_factor = 0.001, clean = clean, verbose = T, mode = "interaction"){

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


  #create binary adjacency matrix first
  dictionary = makeAdjacencyMatrix(tableY = tableY, dict_list = dictionary,
                                   max_factor = max_factor, clean = clean,
                                   verbose = verbose,
                                   mode = mode)

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

#' Starting point to wrangle anansi dictionary list into binary adjacency matrix
#' @description Takes the anansi dictionary in list format and wrangles it into a binary adjacency matrix.
#' @seealso \code{\link{weaveWebFromTables}}
#' @param tableY A matrix containing y features of interest. Rows should be samples and columns should be features. Only used for mode == "interaction".
#' @param dict_list A list that has compound names as entries. For general use, we recommend using the one provided in this package.
#' @param max_factor A proportion (0-1). The maximum proportion of the total dataset that can be a member of a group.
#' @param clean A boolean, whether to clear temporary files during the function.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @param mode A character vector. Can be "interaction" or "membership". Toggles whether to link two datasets based on their interactions or based on shared group membership.
#' @return a binary adjacency matrix with features from tableY as rows and features from tableX as columns.
#'
makeAdjacencyMatrix <- function(dict_list, tableY = NULL, max_factor = 0.001, clean = T, verbose = T, mode = mode){
  if(mode == "interaction"){
    if(verbose){print("Operating in interaction mode")}
    makeAdjacencyMatrixFromList(tableY = tableY, dict_list = dict_list)
  }
  else if(mode == "membership"){
    if(verbose){print("Operating in membership mode")}
    makeAdjacencyMatrixFromGroupMemberList(member_list = dict_list,
                                           max_factor  = max_factor,
                                           clean       = clean,
                                           verbose     = verbose)
  }
}

#' Wrangle anansi dictionary list into binary adjacency matrix
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


#'Wrangle anansi dictionary list into binary adjacency matrix
#' @description Takes the anansi dictionary in list format and wrangles it into a binary adjacency matrix based on group membership.
#' For general use, should probably not be called directly, but rather through \code{\link{weaveWebFromTables}}.
#' @seealso \code{\link{weaveWebFromTables}}
#' @param member_list A list that has compound names as entries. For general use, we recommend using the one provided in this package.
#' @param max_factor A proportion (0-1). The maximum proportion of the total dataset that can be a member of a group.
#' @param clean A boolean, whether to clear temporary files during the function.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a binary adjacency matrix with the group members on both the rows and columns.
#' @importFrom utils combn
#'
makeAdjacencyMatrixFromGroupMemberList <- function(member_list, max_factor = 0.001, clean = T, verbose = T){

  if(verbose){print("Preparing to build adjacency matrix by membership")}

  #Create a list with all pairs of members by group
  memberPairs = lapply(X = member_list, FUN = function(x){combn(x = x, m = 2)})

  #Extract information to create the adjacency matrix
  #This is slow code for large datasets.
  if(verbose){print("Building adjacency matrix...")}

  allMembers <- get_uniques(memberPairs)

  #Create adjacency matrix
  adjMatrix = matrix(nrow = length(allMembers),
                     ncol = length(allMembers),
                     dimnames = list(allMembers, allMembers),
                     data = F)


  if(clean){rm(allMembers); gc()}

  #Populate the empty adjacency matrix with TRUE between features that share a membership
  adjMatrix = fill_mono_adjacency(x = adjMatrix, l = memberPairs,
                                  max_factor = max_factor, clean = clean,
                                  verbose = verbose)

  return(adjMatrix)

}

#' Get unique values from anansi_dict list object
#' @description Helper function to get all unique members in \code{\link{makeAdjacencyMatrixFromGroupMemberList}}
#' @param mp a list of all pairs of members by group.
#' @return a vector of sorted unique member names.
#'
get_uniques <- function(mp){
  ump <- lapply(mp, c)
  ump <- lapply(ump, unique)
  ump <- do.call(c, ump)
  ump <- sort(unique(ump))

  return(ump)
}

#' Populate adjacency matrix based on membership information
#' @description Helper function to populate an adjacency matrix in \code{\link{makeAdjacencyMatrixFromGroupMemberList}}
#' @param x An unfinished adjacency matrix.
#' @param l An element of a list - a matrix with two rows, each row depicting one part of member pairs.
#' @param clean A boolean. Whether to explicitly remove temporary files in order to manage memory.
#' @return A populated adjacency matrix.
#'
assign_adjacency <- function(x, l, clean = T){
  adjmat = x
  adjmat[l[1,],
         l[2,]] <- T
  #And the mirror:
  adjmat[l[2,],
         l[1,]] <- T
  if(clean){
    rm(x)
    gc(verbose = F)}
  return(adjmat)
}


#' Populate adjacency matrix based on membership information. This step can be very heavy on the memory.
#' @description Helper function to populate an adjacency matrix in \code{\link{makeAdjacencyMatrixFromGroupMemberList}}
#' @param x An unfinished adjacency matrix.
#' @param l An element of a list - a matrix with two rows, each row depicting one part of member pairs.
#' @param max_factor A proportion (0-1). The maximum proportion of the total dataset that can be a member of a group.
#' @param clean A boolean. Whether to explicitly remove temporary files in order to manage memory.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return A populated adjacency matrix.
#'
fill_mono_adjacency <- function(x, l, max_factor = 0.001, clean = T, verbose = T){

  members = sapply(l, length) / prod(dim(x)) < max_factor
  stopifnot("No groups have few enough members to satisfy the criteria" = sum(members) >=1)
  if(verbose){
    print(paste(sum(members), "out of", length(members), "groups were exclusive enough to pass the filtering step."))
  }

  l = l[members]

  if(verbose){print("Start populating adjacency matrix.")}

  for(i in 1:length(l)) {
    x = assign_adjacency(x = x, l = l[[i]], clean = clean)
    if(clean){gc()}
    print(i)

  }

  return(x)
}

