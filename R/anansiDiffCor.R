anansiDiffCor = function(web, groups, adjust.method = adjust.method){
  #Extract data from web
  Y              <- web@tableY
  X              <- web@tableX
  PriorKnowledge <- web@dictionary

  #Create result matrices
  out_rvals      <- web@dictionary
  out_pvals      <- !web@dictionary
  out_disjrvals  <- web@dictionary
  out_disjpvals  <- !web@dictionary
  out_emergrvals <- web@dictionary
  out_emergpvals <- !web@dictionary

  #Nested for loop: for every metabolite, do:
  for(y in 1:ncol(Y)){
    #Find the associated functions
    included_functions <- X[,PriorKnowledge[y,] == 1]

    for(x in 1:ncol(included_functions)){
      #For every function associatd to that molecule, fit some models:
      X_tested <- included_functions[,x]

      # fit linear model
      fit      <- lm(Y[,y] ~ X_tested * groups)

      # Calculate p-value for entire model
      fstat    <- summary(fit)$fstatistic
      p        <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

      out_rvals[y, colnames(included_functions)[x]] <- summary(fit)$r.squared
      out_pvals[y, colnames(included_functions)[x]] <- p

      #run ANOVA to determine impact of group on SLOPE of association.
      disj_fit        <- anova(fit)

      #Calculate r.squared
      disj.rsquared  <- disj_fit[3,2] / (disj_fit[3,2] + disj_fit[4,2])

      out_disjrvals[y, colnames(included_functions)[x]]  <- disj.rsquared
      out_disjpvals[y, colnames(included_functions)[x]]  <- disj_fit[3,5]

      #run ANOVA to determine impact of group on STRENGTH of association.
      abs_resid      <- abs(residuals(lm(Y[,y] ~ X_tested)))
      emerg_fit      <- anova(lm(abs_resid ~ groups))

      #Calculate r.squared
      emerg.rsquared <- emerg_fit[1,2] / (emerg_fit[1,2] + emerg_fit[2,2])

      out_emergrvals[y, colnames(included_functions)[x]] <- emerg.rsquared
      out_emergpvals[y, colnames(included_functions)[x]] <- emerg_fit[1,5]

    }
  }
  #Adjust for multiple comparisons
  out_qvals                      <- out_pvals
  out_qvals[PriorKnowledge]      <- p.adjust(out_pvals[PriorKnowledge], method = adjust.method)

  out_disjqvals                  <- out_disjpvals
  out_disjqvals[PriorKnowledge]  <- p.adjust(out_disjpvals[PriorKnowledge], method = adjust.method)

  out_emergqvals                 <- out_emergpvals
  out_emergqvals[PriorKnowledge] <- p.adjust(out_emergpvals[PriorKnowledge], method = adjust.method)

  #Collect into nested list and return results
  return(list(modelfit   = list(r.squared = out_rvals,
                                p         = out_pvals,
                                q         = out_qvals),
              disjointed = list(r.squared = out_disjrvals,
                                p         = out_disjpvals,
                                q         = out_disjqvals),
              emergent   = list(r.squared = out_emergrvals,
                                p         = out_emergpvals,
                                q         = out_emergqvals)))

}

