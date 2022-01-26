anansiCor = function(web, method = "pearson", groups = NULL){
  #Run correlations on subsections of your data
  if(is.null(groups) | any(is.na(groups))){groups = TRUE}
  merge <- cbind(web@tableY[groups,], web@tableX[groups,])
  cors  <- cor(merge, method = method)
  cors_bipartite <- cors[1:ncol(web@tableY), (ncol(web@tableY)+1):ncol(cors)]
  return(cors_bipartite * web@dictionary)
}
