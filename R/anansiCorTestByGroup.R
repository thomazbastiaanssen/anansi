anansiCorTestByGroup = function(web, method = "pearson", groups, adjust.method = adjust.method, verbose = T){
  #Determine all groups
  all_groups = unique(groups)
  
  #Assess correlations for entire data set
  all_out = anansiCorPvalue(web, method = method, adjust.method = adjust.method)
  out_list = list(All = all_out)
  
  #If groups argument is suitable, run subset analysis
  if(length(unique(groups)) < floor(length(groups)/3)){
    
  #If verbose, verbalize. 
  if(verbose){print(paste("Running correlations for the following groups:", 
                          paste(all_groups, collapse = ", "), 
                          "and all together"))}
  
  #Assess correlations for subsets if applicable. 
  for(i in 1:length(all_groups)){
    
    out_list[[i+1]] = anansiCorPvalue(web, method = method, groups = groups == all_groups[i], adjust.method = adjust.method)
  }
  names(out_list)[-1] <- all_groups
  }
  #Return results
  return(out_list)
}
