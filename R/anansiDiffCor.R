#' Run differential correlation analysis for all interacting metabolites and functions.
#' @description Can either take continuous or categorical data for \code{groups}. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param adjust.method Method to adjust p-values for multiple comparisons. \code{adjust.method = "BH"} is the default value. See \code{p.adjust()} in the base R \code{stats} package.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals
#'
anansiDiffCor = function(web, groups, adjust.method = adjust.method, verbose = T){
  #Create a matrix with row and column coordinates to cycle through the relevant comparisons in tableY and tableX.
  which_dictionary <- which(web@dictionary, arr.ind = T, useNames = F)

  stats_out <- apply(which_dictionary, 1, FUN = glm_calc_diff_cor, web = web, groups = groups)

  #Create result container matrices. We take advantage of the fact that true interactions are coded as TRUE, which corresponds to 1,
  #automatically setting all non-canonical interactions as p = 1 and estimate = 0.

  out_rvals      <- web@dictionary
  out_pvals      <- !web@dictionary
  out_disjrvals  <- web@dictionary
  out_disjpvals  <- !web@dictionary
  out_emergrvals <- web@dictionary
  out_emergpvals <- !web@dictionary

  out_rvals[web@dictionary]      <- stats_out[1,]
  out_pvals[web@dictionary]      <- stats_out[2,]
  out_disjrvals[web@dictionary]  <- stats_out[3,]
  out_disjpvals[web@dictionary]  <- stats_out[4,]
  out_emergrvals[web@dictionary] <- stats_out[5,]
  out_emergpvals[web@dictionary] <- stats_out[6,]

  #Adjust for multiple comparisons
  out_qvals                      <- out_pvals
  out_qvals[web@dictionary]      <- p.adjust(out_pvals[web@dictionary], method = adjust.method)

  out_tale        = new("anansiTale",
                        subject    = "model_full",
                        type       = "r.squared",
                        estimates  = out_rvals,
                        p.values   = out_pvals,
                        q.values   = out_qvals)

  out_disjqvals                  <- out_disjpvals
  out_disjqvals[web@dictionary]  <- p.adjust(out_disjpvals[web@dictionary], method = adjust.method)

  out_disjointed = new("anansiTale",
                       subject    = "model_disjointed",
                       type       = "r.squared",
                       estimates  = out_disjrvals,
                       p.values   = out_disjpvals,
                       q.values   = out_disjqvals)

  out_emergqvals                 <- out_emergpvals
  out_emergqvals[web@dictionary] <- p.adjust(out_emergpvals[web@dictionary], method = adjust.method)

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

#' Run differential correlation analysis for all interacting metabolites and functions.
#' @description Should not be run on its own. to be applied by \code{anansiDiffCor()}. Typically, the main \code{anansi()} function will handle this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param which_dictionary A matrix derived from calling \code{which(web@dictionary, arr.ind = T)}. It contains coordinates for the relevant measurements to be compared.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals
#'
glm_calc_diff_cor <- function(web, which_dictionary, groups){
  # Extract relevant values
  y = web@tableY[,which_dictionary[1]]
  x = web@tableX[,which_dictionary[2]]

  vec_out = c(0, 1, 0, 1, 0, 1)

  # fit linear model
  fit      <- lm(y ~ x * groups)

  # Calculate p-value for entire model
  fstat    <- summary(fit)$fstatistic
  p        <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

  vec_out[1] <- summary(fit)$r.squared
  vec_out[2] <- p

  #run ANOVA to determine impact of group on SLOPE of association.
  disj_fit        <- anova(fit)

  #Calculate r.squared
  disj.rsquared  <- disj_fit[3,2] / (disj_fit[3,2] + disj_fit[4,2])

  vec_out[3]  <- disj.rsquared
  vec_out[4]  <- disj_fit[3,5]

  #run ANOVA to determine impact of group on STRENGTH of association.
  abs_resid      <- abs(residuals(lm(y ~ x)))
  emerg_fit      <- anova(lm(abs_resid ~ groups))

  #Calculate r.squared
  emerg.rsquared <- emerg_fit[1,2] / (emerg_fit[1,2] + emerg_fit[2,2])

  vec_out[5] <- emerg.rsquared
  vec_out[6] <- emerg_fit[1,5]

  return(vec_out)
}
