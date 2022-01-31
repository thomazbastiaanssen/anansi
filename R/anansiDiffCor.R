#' Run differential correlation analysis for all interacting metabolites and functions.
#' @description Can either take continuous or categorical data for \code{groups}. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param adjust.method Method to adjust p-values for multiple comparisons. \code{adjust.method = "BH"} is the default value. See \code{p.adjust()} in the base R \code{stats} package.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals
#'
anansiDiffCor = function(web, groups, adjust.method = adjust.method){
  #Extract data from web
  Y              <- web@tableY
  X              <- web@tableX
  PriorKnowledge <- web@dictionary

  #Create result container matrices. We take advantage of the fact that true interactions are coded as TRUE, which corresponds to 1,
  #automatically setting all non-canonical interactions as p = 1 and estimate = 0.

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
      #For every function associated to that molecule, fit some models:
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

  out_tale        = new("anansiTale",
                       subject    = "model_full",
                       type       = "r.squared",
                       estimates  = out_rvals,
                       p.values   = out_pvals,
                       q.values   = out_qvals)

  out_disjqvals                  <- out_disjpvals
  out_disjqvals[PriorKnowledge]  <- p.adjust(out_disjpvals[PriorKnowledge], method = adjust.method)

  out_disjointed = new("anansiTale",
                       subject    = "model_disjointed",
                       type       = "r.squared",
                       estimates  = out_disjrvals,
                       p.values   = out_disjpvals,
                       q.values   = out_disjqvals)

  out_emergqvals                 <- out_emergpvals
  out_emergqvals[PriorKnowledge] <- p.adjust(out_emergpvals[PriorKnowledge], method = adjust.method)

  out_emergent   = new("anansiTale",
                       subject    = "model_emergent",
                       type       = "r.squared",
                       estimates  = out_emergrvals,
                       p.values   = out_emergpvals,
                       q.values   = out_emergqvals)

  #Collect into nested list and return results
  return(list(modelfit   = out_tale,
              disjointed = out_disjointed,
              emergent   = out_emergent))

}

