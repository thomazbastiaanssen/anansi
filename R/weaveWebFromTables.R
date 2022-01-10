weaveWebFromTables = function(table1, table2, dictionary, verbose = T){
  #For conventional use, table 1 should be metabolites and table 2 funcitons. 
  available_table1 = sort(intersect(colnames(table1), rownames(dictionary)))
  available_table2 = sort(intersect(colnames(table2), colnames(dictionary)))
  
  if(verbose){
    print(paste(length(available_table1), "were matched between table 1 and the columns of the adjacency matrix"))
    print(paste(length(available_table2), "were matched between table 2 and the rows of the adjacency matrix"))
  }
  #Select the relevant part of the adjacency matrix
  dictionary = dictionary[available_table1, available_table2]
  
  #Ensure there are no metabolites or functions that never interact
  dictionary = dictionary[rowSums(dictionary) > 0,colSums(dictionary) > 0]
  
  #use the dictionary to clean input tables
  table1 = table1[,rownames(dictionary)]
  table2 = table2[,colnames(dictionary)]
  
  #Return a list with three objects: metabolites, functions and adjacency matrix
  return(list(table1     = table1, 
              table2     = table2, 
              dictionary = dictionary))
}
