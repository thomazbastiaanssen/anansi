#' Run differential correlation analysis for all interacting metabolites and functions.
#' @description Can either take continuous or categorical data for \code{groups}. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param formula A formula object. Used to assess differential associations.
#' @param reff A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.
#' @param modeltype A string, either "lm" or "lmer" depending on the type of model that should be ran.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals terms
#' @importFrom future.apply future_apply
#'
anansiDiffCor = function(web, metadata, groups, formula, reff, modeltype, verbose = T){
  #Create a matrix with row and column coordinates to cycle through the relevant comparisons in tableY and tableX.
  which_dictionary <- which(web@dictionary, arr.ind = T, useNames = F)

  all_terms <- attr(terms(formula), "term.labels")

  stats_out <- future.apply::future_apply(which_dictionary, 1, FUN = model_picker,
                                          web = web, formula = formula, metadata = metadata,
                                          reff = reff, modeltype = modeltype)

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

  #Adjust for multiple comparisons
  out_qvals                      <- out_pvals
  out_qvals[web@dictionary]      <- NA

  out_tale  <- list(full = new("anansiTale",
                               subject    = "model_full",
                               type       = "r.squared",
                               estimates  = out_rvals,
                               p.values   = out_pvals,
                               q.values   = out_qvals))

  out_disjointed = vector(mode =  "list", length = length(all_terms))

  for(t in 1:length(all_terms)){
    out_disjrvals[web@dictionary]  <- stats_out[1 + t * 2,]
    out_disjpvals[web@dictionary]  <- stats_out[2 + t * 2,]

    #Adjust for multiple comparisons
    out_disjqvals                  <- out_disjpvals
    out_disjqvals[web@dictionary]  <- NA


    out_disjointed[[t]] <- new("anansiTale",
                          subject    = paste("model_disjointed", all_terms[t], sep = "_"),
                          type       = "r.squared",
                          estimates  = out_disjrvals,
                          p.values   = out_disjpvals,
                          q.values   = out_disjqvals)
  }
  names(out_disjointed) <- all_terms

  out_emergent = vector(mode =  "list", length = length(all_terms))

  for(t in 1:length(all_terms)){
    out_emergrvals[web@dictionary]  <- stats_out[1 + 2 * length(all_terms) + t * 2,]
    out_emergpvals[web@dictionary]  <- stats_out[2 + 2 * length(all_terms) + t * 2,]

    #Adjust for multiple comparisons
    out_emergqvals                  <- out_emergpvals
    out_emergqvals[web@dictionary]  <- NA


    out_emergent[[t]] <- new("anansiTale",
                             subject    = paste("model_emergent", all_terms[t], sep = "_"),
                             type       = "r.squared",
                             estimates  = out_emergrvals,
                             p.values   = out_emergpvals,
                             q.values   = out_emergqvals)
  }
  names(out_emergent) <- all_terms


  #Collect into nested list and return results
  return(list(modelfit   = out_tale,
              disjointed = out_disjointed,
              emergent   = out_emergent))

}

#' Choose the appropriate model call for \code{anansiDiffCor()}.
#' @description Should not be run on its own. to be applied by \code{anansiDiffCor()}. Typically, the main \code{anansi()} function will handle this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param which_dictionary A matrix derived from calling \code{which(web@dictionary, arr.ind = T)}. It contains coordinates for the relevant measurements to be compared.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param formula A formula object. Used to assess differential associations.
#' @param reff A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.
#' @param modeltype A string, either "lm" or "lmer" depending on the type of model that should be ran.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#'
model_picker <- function(web, which_dictionary, metadata, formula, reff = NULL, modeltype = "lm"){
  if(modeltype == "lm"){
    if(identical(web@tableY, web@tableX)){
      return(glm_calc_diff_prop(web = web, which_dictionary = which_dictionary, metadata = metadata))
    }
    if(class(web) == "argonansiWeb"){
      return(glm_argonaut_calc_diff_cor(web = web, which_dictionary = which_dictionary, formula = formula, metadata = metadata))
    }

    return(glm_calc_diff_cor(web = web, which_dictionary = which_dictionary, formula = formula, metadata = metadata))
  }
  if(modeltype == "lmer"){
    if(identical(web@tableY, web@tableX)){
      return(glmer_calc_diff_prop(web = web, which_dictionary = which_dictionary, metadata = metadata, reff = reff))
    }
    if(class(web) == "argonansiWeb"){
      return(glmer_argonaut_calc_diff_cor(web = web, which_dictionary = which_dictionary, formula = formula, metadata = metadata, reff = reff))
    }

    return(glmer_calc_diff_cor(web = web, which_dictionary = which_dictionary, formula = formula, metadata = metadata, reff = reff))
  }
}

#' Run differential correlation analysis for all interacting metabolites and functions.
#' @description Should not be run on its own. to be applied by \code{anansiDiffCor()}. Typically, the main \code{anansi()} function will handle this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param which_dictionary A matrix derived from calling \code{which(web@dictionary, arr.ind = T)}. It contains coordinates for the relevant measurements to be compared.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param formula A formula object. Used to assess differential associations.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals na.exclude
#'
glm_calc_diff_cor <- function(web, which_dictionary, metadata, formula){

  meta_glm = metadata

  # Extract relevant values
  meta_glm$y = web@tableY[,which_dictionary[1]]
  meta_glm$x = web@tableX[,which_dictionary[2]]

  all_terms <- attr(terms(formula), "term.labels")

  vec_out = rep(c(0, 1), 1 + 2 * length(all_terms))
  formula_full = update.formula(formula, y ~ x * 1 * (.))

  # fit linear model
  fit      <- lm(formula = formula_full, data = meta_glm)

    # Calculate p-value for entire model
  fstat    <- summary(fit)$fstatistic
  p        <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

  vec_out[1] <- summary(fit)$r.squared
  vec_out[2] <- p

  #run ANOVA to determine impact of group on SLOPE of association.
  disj_fit        <- anova(fit)

  #Calculate r.squared
  #get the residual sum of squares
  resid.ss                 <- sum(residuals(fit)^2)
  target_disj_interactions <- grep("^x:",row.names(disj_fit))
  disj.rsquared            <- disj_fit[target_disj_interactions,2] / (disj_fit[target_disj_interactions,2] + resid.ss)

  vec_out[1 + (1:length(all_terms))*2]  <- disj.rsquared
  vec_out[2 + (1:length(all_terms))*2]  <- disj_fit[target_disj_interactions,5]

  #run ANOVA to determine impact of group on STRENGTH of association.
  formula_emerg = update.formula(abs_resid ~ ., formula)

  meta_glm$abs_resid      <- abs(residuals(lm(y ~ x, data = meta_glm, na.action = na.exclude)))

  emerg_fit      <- lm(formula = formula_emerg, data = meta_glm)
  emerg_anova    <- anova(emerg_fit)

  #Calculate r.squared
  resid_emerg.ss <- sum(residuals(emerg_fit)^2)
  target_emerg_interactions <- 1:length(all_terms)

  emerg.rsquared <- emerg_anova[target_emerg_interactions,2] / (emerg_anova[target_emerg_interactions,2] + resid_emerg.ss)

  vec_out[1 + 2 * length(all_terms) + (1:length(all_terms))*2] <- emerg.rsquared
  vec_out[2 + 2 * length(all_terms) + (1:length(all_terms))*2] <- emerg_anova[target_emerg_interactions,5]

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
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param formula A formula object. Used to assess differential associations.
#' @param reff A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals fitted cor na.exclude as.formula
#' @importFrom lme4 lmer
#
glmer_calc_diff_cor <- function(web, which_dictionary, metadata, formula, reff){

  meta_glmer = metadata

  # Extract relevant values
  meta_glmer$y = web@tableY[,which_dictionary[1]]
  meta_glmer$x = web@tableX[,which_dictionary[2]]

  formula_full = update.formula(formula, y ~ x * 1 * (.))

  which_terms <- grep(pattern = "^x:", attr(terms(formula_full), "term.labels"))
  all_terms <- attr(terms(formula_full), "term.labels")[which_terms]
  #param      r, p, r, p, r, p
  vec_out = rep(c(0, 1), 1 + 2 * length(all_terms))

  # fit null model to compute p-values later
  f.null   <- lm(y ~ 1, data = meta_glmer, na.action = na.exclude)
  # fit complete mixed effect model
  formula_full = update.formula(formula_full, .~. + (1|reff))

  f.full   <- suppressMessages(lme4::lmer(formula = formula_full, data = meta_glmer, REML = F, na.action = na.exclude))

  #Compare null to full model for get a p-value
  f_comp   <- anova.merMod(f.null, f.full)
  p        <- f_comp["f.full", "Pr(>Chisq)"]

  #Here, we calculate R^2 for the complete fitted model using pearson's method
  vec_out[1] <- cor(meta_glmer$y, fitted(f.full), method = "pearson", use = "pairwise.complete.obs")^2
  vec_out[2] <- p

  target_disj_interactions <- attr(terms(formula_full), "term.labels")[grep(pattern = "^x:", attr(terms(formula_full), "term.labels"))]

  #Fit a separate null model to compute effect of groups value on slope (ie interaction).
  target.list <- lapply(target_disj_interactions, FUN = function(x){as.formula(paste0(".~.-", x))})

  #update to remove each term sequentially
  f.null.list <- lapply(target.list, FUN = function(x){update.formula(formula_full, new =x )})

  #fit models
  f.null.list <- lapply(X = f.null.list, FUN = function(f){suppressMessages(lme4::lmer(as.formula(f),data = meta_glmer, REML = F, na.action = na.exclude))})

  names(f.null.list) <- target_disj_interactions
  #compare null to full model for get a p-value

  f_comp2   <- lapply(f.null.list, function(x){anova.merMod(x, f.full)["f.full", "Pr(>Chisq)"]})

  f_comp2.rsq.list <- lapply(f.null.list, function(f){cor(residuals(f), fitted(f.full), method = "pearson", use = "pairwise.complete.obs")^2 })
  #Here, we calculate R^2 for the complete fitted model using pearson's method
  vec_out[1 + (1:length(all_terms))*2] <- unlist(f_comp2.rsq.list)
  vec_out[2 + (1:length(all_terms))*2] <- unlist(f_comp2)


  #run ANOVA to determine impact of group on STRENGTH of association.
  meta_glmer$abs_resid      <- abs(residuals(suppressMessages(lme4::lmer(y ~ x + (1|reff), data = meta_glmer, REML = F, na.action = na.exclude))))




  #run ANOVA to determine impact of group on STRENGTH of association.
  formula_emerg = update.formula(abs_resid ~ ., formula)

  #exactly the same as the lm version for here on
  emerg_fit      <- lm(formula = formula_emerg, data = meta_glmer)
  emerg_anova    <- anova(emerg_fit)

  #Calculate r.squared
  resid_emerg.ss <- sum(residuals(emerg_fit)^2)
  target_emerg_interactions <- 1:length(all_terms)

  emerg.rsquared <- emerg_anova[target_emerg_interactions,2] / (emerg_anova[target_emerg_interactions,2] + resid_emerg.ss)

  vec_out[1 + 2 * length(all_terms) + (1:length(all_terms))*2] <- emerg.rsquared
  vec_out[2 + 2 * length(all_terms) + (1:length(all_terms))*2] <- emerg_anova[target_emerg_interactions,5]

  return(vec_out)
}


#' Run differential association analysis within a dataset.
#' @description Should not be run on its own. to be applied by \code{anansiDiffCor()}. Typically, the main \code{anansi()} function will handle this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param which_dictionary A matrix derived from calling \code{which(web@dictionary, arr.ind = T)}. It contains coordinates for the relevant measurements to be compared.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param formula A formula object. Used to assess differential associations.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals na.exclude
#'
glm_calc_diff_prop <- function(web, which_dictionary, metadata, formula){

  meta_glm = metadata
  all_terms <- attr(terms(formula), "term.labels")

  # Extract relevant values
  y = web@tableY[,which_dictionary[1]]
  x = web@tableX[,which_dictionary[2]]

  meta_glm$prop = log(y/x)

  vec_out = rep(c(0, 1), 1 + 2 * length(all_terms))

  # fit linear model
  fit      <- lm(update.formula(formula, prop ~ .), data = meta_glm)

  # Calculate p-value for entire model
  fstat    <- summary(fit)$fstatistic
  p        <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

  vec_out[1] <- summary(fit)$r.squared
  vec_out[2] <- p

  #run ANOVA to determine impact of group on SLOPE of association.
  disj_fit        <- anova(fit)

  #Calculate r.squared
  #get the residual sum of squares
  resid.ss                 <- sum(residuals(fit)^2)
  target_disj_interactions <- 1:length(all_terms)
  disj.rsquared            <- disj_fit[target_disj_interactions,2] / (disj_fit[target_disj_interactions,2] + resid.ss)

  vec_out[1 + (1:length(all_terms))*2]  <- disj.rsquared
  vec_out[2 + (1:length(all_terms))*2]  <- disj_fit[target_disj_interactions,5]


  disj.rsquared  <- disj_fit[1,2] / (disj_fit[1,2] + disj_fit[2,2])


  #run ANOVA to determine impact of group on STRENGTH of association.
  formula_emerg = update.formula(abs_resid ~ ., formula)

  meta_glm$abs_resid  <- abs(residuals(lm(prop ~ 1, na.action = na.exclude, data = meta_glm)))

  emerg_fit      <- lm(formula = formula_emerg, data = meta_glm)
  emerg_anova    <- anova(emerg_fit)

  #Calculate r.squared
  resid_emerg.ss <- sum(residuals(emerg_fit)^2)
  target_emerg_interactions <- 1:length(all_terms)

  emerg.rsquared <- emerg_anova[target_emerg_interactions,2] / (emerg_anova[target_emerg_interactions,2] + resid_emerg.ss)

  vec_out[1 + 2 * length(all_terms) + (1:length(all_terms))*2] <- emerg.rsquared
  vec_out[2 + 2 * length(all_terms) + (1:length(all_terms))*2] <- emerg_anova[target_emerg_interactions,5]

  return(vec_out)
}

#' Run differential association analysis within a dataset.
#' @description Should not be run on its own. to be applied by \code{anansiDiffCor()}. Typically, the main \code{anansi()} function will handle this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param which_dictionary A matrix derived from calling \code{which(web@dictionary, arr.ind = T)}. It contains coordinates for the relevant measurements to be compared.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param formula A formula object. Used to assess differential associations.
#' @param reff A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals fitted cor na.exclude as.formula
#' @importFrom lme4 lmer
#
glmer_calc_diff_prop <- function(web, which_dictionary, metadata, formula, reff){

  meta_glmer = metadata
  all_terms <- attr(terms(formula), "term.labels")

  # Extract relevant values
  y = web@tableY[,which_dictionary[1]]
  x = web@tableX[,which_dictionary[2]]

  prop = log(y/x)
  #param      r, p, r, p, r, p
  vec_out = rep(c(0, 1), 1 + 2 * length(all_terms))

  # fit null model to compute p-values later
  f.null   <- lm(prop ~ 1, na.action = na.exclude, data = meta_glmer)

  # fit complete mixed effect model
  formula_full = update.formula(formula, prop ~ . )

  formula_full = update.formula(formula_full, .~ . + (1|reff))

  f.full   <- suppressMessages(lme4::lmer(formula = formula_full, data = meta_glmer, REML = F, na.action = na.exclude))

  #Compare null to full model for get a p-value
  f_comp   <- anova.merMod(f.null, f.full)
  p        <- f_comp["f.full", "Pr(>Chisq)"]

  #Here, we calculate R^2 for the complete fitted model using pearson's method
  vec_out[1] <- cor(prop, fitted(f.full), method = "pearson", use = "pairwise.complete.obs")^2
  vec_out[2] <- p


  target_disj_interactions <- all_terms

  #Fit a separate null model to compute effect of groups value on slope (ie interaction).
  target.list <- lapply(target_disj_interactions, FUN = function(x){as.formula(paste0(".~.-", x ))})

  #update to remove each term sequentially
  f.null.list <- lapply(target.list, FUN = function(x){update.formula(formula_full, new =x )})

  #fit models
  f.null.list <- lapply(X = f.null.list, FUN = function(f){suppressMessages(lme4::lmer(as.formula(f),data = meta_glmer, REML = F, na.action = na.exclude))})

  names(f.null.list) <- target_disj_interactions
  #compare null to full model for get a p-value

  f_comp2   <- lapply(f.null.list, function(x){anova.merMod(x, f.full)["f.full", "Pr(>Chisq)"]})

  f_comp2.rsq.list <- lapply(f.null.list, function(f){cor(residuals(f), fitted(f.full), method = "pearson", use = "pairwise.complete.obs")^2 })
  #Here, we calculate R^2 for the complete fitted model using pearson's method
  vec_out[1 + (1:length(all_terms))*2] <- unlist(f_comp2.rsq.list)
  vec_out[2 + (1:length(all_terms))*2] <- unlist(f_comp2)


  #run ANOVA to determine impact of group on STRENGTH of association.
  meta_glmer$abs_resid      <- abs(residuals(suppressMessages(lme4::lmer(prop ~ 1 + (1|reff), data = meta_glmer, REML = F, na.action = na.exclude))))

  #run ANOVA to determine impact of group on STRENGTH of association.
  formula_emerg = update.formula(abs_resid ~ ., formula)

  #exactly the same as the lm version for here on
  emerg_fit      <- lm(formula = formula_emerg, data = meta_glmer)
  emerg_anova    <- anova(emerg_fit)

  #Calculate r.squared
  resid_emerg.ss <- sum(residuals(emerg_fit)^2)
  target_emerg_interactions <- 1:length(all_terms)

  emerg.rsquared <- emerg_anova[target_emerg_interactions,2] / (emerg_anova[target_emerg_interactions,2] + resid_emerg.ss)

  vec_out[1 + 2 * length(all_terms) + (1:length(all_terms))*2] <- emerg.rsquared
  vec_out[2 + 2 * length(all_terms) + (1:length(all_terms))*2] <- emerg_anova[target_emerg_interactions,5]

  return(vec_out)
}
