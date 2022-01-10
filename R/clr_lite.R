impute_zeroes = function(vec, method = "logunif"){
  if(! method %in% c("logunif", "unif", "const")){stop("`method` must be exactly `logunif`, `unif` or `const`")}
  
  #Find detection limit
  DL = min(vec[vec != 0])
  if(method == "logunif"){
    vec[vec == 0] = DL/(10^(runif(n = sum(vec == 0), min =  0, max = 1)))
  }
  else if(method == "unif"){
    vec[vec == 0] = runif(n = sum(vec == 0), min =  0.1*DL, max = DL)
  }
  else if(method == "const"){
    vec[vec == 0] = 0.65 * DL
  }
  return(vec)
} 


clr_imputed = function(vec, method = "logunif", replicates = 1000){
  if(! method %in% c("logunif", "unif", "const")){stop("`method` must be exactly `logunif`, `unif` or `const`")}
  return(apply(replicate(replicates, compositions::clr(impute_zeroes(vec = vec, method = method))), 1, median))
}


clr_lite = function(counts, samples_are = "cols", method = "logunif", replicates = 1000) 
{
  temp_counts = counts
  
  if(! method %in% c("logunif", "unif", "const"))
  {stop("`method` must be exactly `logunif`, `unif` or `const`")}
  
  if(samples_are == "rows"){
    temp_counts = data.frame(t(temp_counts))
  }
  
  temp_counts = apply(X          = temp_counts, 
                      MARGIN     = 2, 
                      FUN        = clr_imputed, 
                      method     = method, 
                      replicates = replicates)
  
  if(samples_are == "rows"){
    temp_counts = data.frame(t(temp_counts))
  }
  
  clr_counts = data.frame(temp_counts)
  rownames(clr_counts) = rownames(counts)
  colnames(clr_counts) = colnames(counts)
  return(clr_counts)
}
