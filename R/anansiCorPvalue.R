anansiCorPvalue = function(web, method = "pearson", groups = NULL, adjust.method = adjust.method) {
  #Compute correlation coefficients
  r    <- anansiCor(web = web, method = method, groups = groups)
  
  if(is.null(groups)){groups = TRUE}
  #Compute t-statistics based on the n and the correlation coefficient
  n    <- web$dictionary
  n[T] <- nrow(web$table1[groups,])
  t    <- (r*sqrt(n-2))/sqrt(1-r^2)
  
  #Compute p-values based on t and n. 
  p    <- 2*(1 - pt(abs(t),(n-2)))
  
  #Compute naive adjusted p-values 
  q    <- p
  q[web$dictionary] <- p.adjust(p[web$dictionary], method = adjust.method)
  
  #Collate correlation coefficients, p-values and q-values into list
  out  <- list(r, p, q)
  names(out) <- c("r", "p", "q")
  return(out)
}
