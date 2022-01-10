anansi = function(web, method = "pearson", groups, adjust.method = "BH", verbose = T){
  
  if(verbose){print("Running annotation-based correlations")}
  out_cors  = anansiCorTestByGroup(web = web, method = method, groups = groups, adjust.method = adjust.method, verbose = verbose)
  
  if(verbose){print("Fitting models for differential correlation testing")}
  out_models = anansiDiffCor(web = web, groups = groups, adjust.method = adjust.method)

  return(list(cor_results = out_cors, model_results = out_models))
}
