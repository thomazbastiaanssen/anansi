#' Impute zeroes and perform a centered log-ratio (CLR) transformation on your count table
#' @param counts A compositional count table.
#' @param samples_are Either "cols" or "rows". Default is "cols". Denotes whether the columns or rows depict indicidual samples.
#' @param method The method for zero imputation. One of "logunif", "unif" or "const".
#' @param replicates An integer. for the two random sampling methods, if this is larger than 1, every zero will be imputed that many times. The median of the CLR of all those replicates will be returned.
#' @return A CLR-transformed count table.
#' @export
#'
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

#' compute CLR using Aitchison's method
#' @param x A vector of compositional data without zeroes.
#' @return a vector of CLR-transformed data
#'
anansi_compute_clr = function(x){
  #compute CLR using Aitchison's method
  return(log(x/exp(mean(log(x)))))
}

#' Replace zeroes with non-zero values in order to perform a CLR-transformation
#' @param vec A vector of compositional data that may contain zeroes.
#' @param method The method for zero imputation. One of "logunif", "unif" or "const".
#' @return a vector with all the zeroes replaced with non-zero values.
#'
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

#' Resample random values, perform CLR over each iteration and return the median result.
#' @param vec A vector of compositional data that may contain zeroes.
#' @param method The method for zero imputation. One of "logunif", "unif" or "const".
#' @param replicates A positive integer. Default is 1000. Controls how many replicates the median should be taken over.
#' @return a vector of CLR-transformed data
#'
clr_imputed = function(vec, method = "logunif", replicates = 1000){
  if(! method %in% c("logunif", "unif", "const")){stop("`method` must be exactly `logunif`, `unif` or `const`")}
  return(apply(replicate(replicates, anansi_compute_clr(impute_zeroes(vec = vec, method = method))), 1, median))
}



