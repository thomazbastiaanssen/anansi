anansiCorTestByGroup = function(web, method = "pearson", groups, adjust.method = adjust.method, verbose = T){
  #Determine all groups
  all_groups = unique(groups)

  #Assess correlations for entire data set
  all_out = anansiCorPvalue(web, method = method, adjust.method = adjust.method)
  all_out@type = paste(all_out@type, "All", sep = "_")
  out_list = list(All = all_out)

  #If groups argument is suitable, run subset analysis
  if(length(unique(groups)) < floor(length(groups)/3)){

  #If verbose, verbalize.
  if(verbose){print(paste("Running correlations for the following groups:",
                          paste(all_groups, collapse = ", "),
                          "and all together"))}

  #Assess correlations for subsets if applicable.
  for(i in 1:length(all_groups)){
    out_by_group = anansiCorPvalue(web, method = method, groups = groups == all_groups[i], adjust.method = adjust.method)
    out_by_group@type = paste(out_by_group@type, all_groups[i], sep = "_")

    out_list[[i+1]] = out_by_group
  }
  names(out_list)[-1] <- all_groups
  }
  #Return results
  return(out_list)
}
