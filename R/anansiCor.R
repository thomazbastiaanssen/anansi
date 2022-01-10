anansiCor = function(web, method = "pearson", groups = NULL){
  #Run correlations on subsections of your data
  if(is.null(groups)){groups = TRUE}
  merge <- cbind(web$table1[groups,], web$table2[groups,])
  cors  <- cor(merge, method = method)
  cors_bipartite <- cors[1:ncol(web$table1), (ncol(web$table1)+1):ncol(cors)]
  return(cors_bipartite * web$dictionary)
}
