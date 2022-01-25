weaveWebFromTables = function(tableY, tableX, dictionary, verbose = T){
  #For conventional use, table Y should be metabolites and table X functions.

  # The two 'table' data.frames MUST have row
  # and column names that are unique, and
  # look like the following:
  #
  #               f1  f2  f3  f4  f5  f6
  #   sample1      0   0   2   0   0   1
  #   sample2     20   8  12   5  19  26
  #   sample3      3   0   2   0   0   0
  #       ... many more rows ...
  #
  # the column names need to correspond to the names in the dictionary.

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

  #Return an anansiWeb with three objects: typically metabolites, functions and adjacency matrix
  return(new("anansiWeb",
             tableY     = tableY,
             tableX     = tableX,
             dictionary = dictionary))
}
