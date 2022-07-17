#' Run differential correlation analysis for all interacting metabolites and functions.
#' @description Can either take continuous or categorical data for \code{groups}. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param reff A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.
#' @param modeltype A string, either "lm" or "lmer" depending on the type of model that should be ran.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals
#' @importFrom future.apply future_apply
#'
anansiDiffCor = function(web, groups, reff, modeltype, verbose = T){
  #Create a matrix with row and column coordinates to cycle through the relevant comparisons in tableY and tableX.
  which_dictionary <- which(web@dictionary, arr.ind = T, useNames = F)

  stats_out <- future.apply::future_apply(which_dictionary, 1, FUN = model_picker,
                                          web = web, groups = groups, reff = reff, modeltype = modeltype)

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
  out_qvals[web@dictionary]      <- NA

  out_tale        = new("anansiTale",
                        subject    = "model_full",
                        type       = "r.squared",
                        estimates  = out_rvals,
                        p.values   = out_pvals,
                        q.values   = out_qvals)

  out_disjqvals                  <- out_disjpvals
  out_disjqvals[web@dictionary]  <- NA

  out_disjointed = new("anansiTale",
                       subject    = "model_disjointed",
                       type       = "r.squared",
                       estimates  = out_disjrvals,
                       p.values   = out_disjpvals,
                       q.values   = out_disjqvals)

  out_emergqvals                 <- out_emergpvals
  out_emergqvals[web@dictionary] <- NA

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

#' Choose the appropriate model call for \code{anansiDiffCor()}.
#' @description Should not be run on its own. to be applied by \code{anansiDiffCor()}. Typically, the main \code{anansi()} function will handle this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param which_dictionary A matrix derived from calling \code{which(web@dictionary, arr.ind = T)}. It contains coordinates for the relevant measurements to be compared.
#' @param groups A categorical or continuous vector necessary for differential correlations. Typically a state or treatment score.
#' @param reff A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.
#' @param modeltype A string, either "lm" or "lmer" depending on the type of model that should be ran.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#'
model_picker <- function(web, which_dictionary, groups, reff = NULL, modeltype = "lm"){
  if(modeltype == "lm"){
    return(glm_calc_diff_cor(web = web, which_dictionary = which_dictionary, groups = groups))
  }
  if(modeltype == "lmer"){
    return(glmer_calc_diff_cor(web = web, which_dictionary = which_dictionary, groups = groups, reff = reff))
  }
}

#' Run differential correlation analysis for all interacting metabolites and functions.
#' @description Should not be run on its own. to be applied by \code{anansiDiffCor()}. Typically, the main \code{anansi()} function will handle this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param which_dictionary A matrix derived from calling \code{which(web@dictionary, arr.ind = T)}. It contains coordinates for the relevant measurements to be compared.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals na.exclude
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
  abs_resid      <- abs(residuals(lm(y ~ x, na.action = na.exclude)))
  emerg_fit      <- anova(lm(abs_resid ~ groups))

  #Calculate r.squared
  emerg.rsquared <- emerg_fit[1,2] / (emerg_fit[1,2] + emerg_fit[2,2])

  vec_out[5] <- emerg.rsquared
  vec_out[6] <- emerg_fit[1,5]

  return(vec_out)
}

#' Import anova() for merMod objects from lme4 to compare lm vs lmer with null model first.
#' @description The lme4:::anova.merMod method for anova.
#' @param object merMod object
#' @param ...   further such objects
#' @param refit should objects be refitted with ML (if applicable)
#' @param model.names character vectors of model names to be used in the anova table.
#' @return an "anova" data frame; the traditional (S3) result of anova()
#' @importFrom utils getFromNamespace
#' @importFrom stats anova
#
anova.merMod <- utils::getFromNamespace("anova.merMod", "lme4")


#' Run differential correlation analysis for all interacting metabolites and functions.
#' @description Should not be run on its own. to be applied by \code{anansiDiffCor()}. Typically, the main \code{anansi()} function will handle this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param which_dictionary A matrix derived from calling \code{which(web@dictionary, arr.ind = T)}. It contains coordinates for the relevant measurements to be compared.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param reff A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals fitted cor na.exclude
#' @importFrom lme4 lmer
#
glmer_calc_diff_cor <- function(web, which_dictionary, groups, reff){
  # Extract relevant values
  y = web@tableY[,which_dictionary[1]]
  x = web@tableX[,which_dictionary[2]]

  #param      r, p, r, p, r, p
  vec_out = c(0, 1, 0, 1, 0, 1)

  # fit null model to compute p-values later
  f.null   <- lm(y ~ 1, na.action = na.exclude)
  # fit complete mixed effect model
  f.full   <- suppressMessages(lme4::lmer(y ~ x * groups + (1|reff), REML = F, na.action = na.exclude))

  #Compare null to full model for get a p-value
  f_comp   <- anova.merMod(f.null, f.full)
  p        <- f_comp["f.full", "Pr(>Chisq)"]

  #Here, we calculate R^2 for the complete fitted model using pearson's method
  vec_out[1] <- cor(y, fitted(f.full),method = "pearson", use = "pairwise.complete.obs")^2
  vec_out[2] <- p

  #Fit a separate null model to compute effect of groups value on slope (ie interaction).
  f.null2  <- suppressMessages(lme4::lmer(y ~ x + groups + (1|reff), REML = F, na.action = na.exclude))

  #compare null to full model for get a p-value
  f_comp2   <- anova.merMod(f.null2, f.full)
  p2        <- f_comp2["f.full", "Pr(>Chisq)"]

  #Here, we calculate R^2 for the complete fitted model using pearson's method
  vec_out[3] <- cor(residuals(f.null2), fitted(f.full),method = "pearson", use = "pairwise.complete.obs")^2
  vec_out[4] <- p2


  #run ANOVA to determine impact of group on STRENGTH of association.
  abs_resid      <- abs(residuals(suppressMessages(lme4::lmer(y ~ x + (1|reff), REML = F, na.action = na.exclude))))

  #exactly the same as the lm version for here on
  emerg_fit      <- anova(lm(abs_resid ~ groups))

  #Calculate r.squared
  emerg.rsquared <- emerg_fit[1,2] / (emerg_fit[1,2] + emerg_fit[2,2])

  vec_out[5] <- emerg.rsquared
  vec_out[6] <- emerg_fit[1,5]

  return(vec_out)
}
