#' An S4 class to contain all metabolomics and functional input data as well as a dictionary to link them.
#' Can only be used with the argonaut library.
#'
#' @slot tableY A matrix of metabolomics data. Rows are samples and columns are features.
#' @slot tableX A stratified feature table containing features of interest. dim1 should be samples, dim2 should be features and dim3 should be subtypes.
#' @slot dictionary A binary adjacency matrix. Typically generated using the \code{weaveWebFromTables()} function.
#' @description argonansiWeb is the main container that will hold your input data thoughout the \code{anansi} pipeline.
#' @importClassesFrom  argonaut stratifiedFeatureTable
#'
setClass("argonansiWeb",
  slots = c(
    tableX.sft   = "stratifiedFeatureTable",
    strat_dict   = "matrix"
  ), contains = "anansiWeb"
)



#' Create an argonansiWeb object from two 'omics tables and a dictionary
#' @description This function will take two tables of 'omics data, for example metabolomics and stratified functional microbiome data. It will also take a dictionary list, which is provided in this package.
#' @param tableY A table containing features of interest. Rows should be samples and columns should be features. The Y and X refer to the position of the features in a formula: Y ~ X.
#' @param stratifiedTableX A stratified feature table containing features of interest. dim1 should be samples, dim2 should be features and dim3 should be subtypes.
#' @param dictionary A list that has feature names from tableY as names. Default is the dictionary provided in this package.
#' @param mode A character vector. Can be "interaction" or "membership". Toggles whether to link two datasets based on their interactions or based on shared group membership.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @param prune A boolean. Toggles whether to prune particularly large groups.
#' @param max_sds A numeric. Only relevant for prune == TRUE. How many SDs larger than the median a group can be before it is pruned.
#' For general use, we recommend sticking to that one. You can access the dictionary like this: \code{data(dictionary)}
#' @return an \code{anansiWeb} object. Web is used as input for most of the main workflow of anansi.
#' @export
weaveWebFromStratifiedTables <- function(tableY, stratifiedTableX, dictionary = anansi::anansi_dic, verbose = T, mode = "interaction", prune = F, max_sds = 3) {
  stopifnot("stratifiedTableX needs to be a stratifiedFeatureTable object from the argonaut package." = "stratifiedFeatureTable" %in% class(stratifiedTableX))

  tableX <- argonaut::apply_by(X = stratifiedTableX, MARGIN = 3, mean)

  if (length(dictionary) == 1) {
    if (dictionary == "none") {
      if (verbose) {
        print("No dictionary provided, preparing for all vs all analysis. ")
      }
      dictionary <- mock_dictionary(tableY = tableY, tableX = tableX)
    }
  }

  stopifnot("the mode argument needs to be interaction or membership." = mode %in% c("interaction", "membership"))
  # For conventional use, table Y should be metabolites and table X functions.

  # The two 'table' matrices MUST have row
  # and column names that are unique, and
  # look like the following:
  #
  #               f1  f2  f3  f4  f5  f6
  #   sample1      0   0   2   0   0   1
  #   sample2     20   8  12   5  19  26
  #   sample3      3   0   2   0   0   0
  #       ... many more rows ...
  #

  if (prune) {
    dictionary <- prune_dictionary_for_exclusivity(
      dict_list = dictionary,
      max_sds = max_sds, verbose = verbose
    )
  }

  # create binary adjacency matrix first
  dictionary <- makeAdjacencyMatrix(
    tableY = tableY, dict_list = dictionary,
    verbose = verbose, mode = mode
  )

  # If we're looking at single data set, don't do associations with yourself. Set diagonal to FALSE.
  if (identical(tableX, tableY)) {
    diag(dictionary) <- F
  }

  # Check if input tables have the same names and the same length.
  if (!identical(row.names(tableY), row.names(tableX))) {
    warning("The row names of tableY and tableX do not correspond. Please make sure they are in the same order.")
  }
  if (nrow(tableY) != nrow(tableX)) {
    stop("tableY and tableX do not have the same number of rows/observations.")
  }

  available_tableY <- sort(intersect(colnames(tableY), rownames(dictionary)))
  available_tableX <- sort(intersect(colnames(tableX), colnames(dictionary)))

  if (verbose) {
    print(paste(length(available_tableY), "were matched between table 1 and the columns of the adjacency matrix"))
    print(paste(length(available_tableX), "were matched between table 2 and the rows of the adjacency matrix"))
  }
  # Select the relevant part of the adjacency matrix
  dictionary <- dictionary[available_tableY, available_tableX]

  # Ensure there are no features that never interact
  dictionary <- dictionary[rowSums(dictionary) > 0, colSums(dictionary) > 0]

  # use the dictionary to clean input tables
  tableY <- tableY[, rownames(dictionary)]
  # tableX = tableX[,colnames(dictionary)]
  stratifiedTableX <- argonaut::sft(stratifiedTableX[, colnames(dictionary), ])
  tableX <- get_strat_X(stratifiedTableX)
  strat_dict <- as.matrix(dictionary)[, rep(1:ncol(as.matrix(dictionary)), argonaut::apply_by(stratifiedTableX, 3, length)[1, ])]

  # Return an argonansiWeb object with five slots: typically metabolites, functions, stratified functions, a stratified adjacency matrix and aggregated an adjacency matrix
  return(new("argonansiWeb",
    tableY     = as.matrix(tableY),
    tableX     = as.matrix(tableX),
    tableX.sft = stratifiedTableX,
    strat_dict = as.matrix(strat_dict),
    dictionary = as.matrix(dictionary)
  ))
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
glm_argonaut_calc_diff_cor <- function(web, which_dictionary, metadata, formula) {
  meta_glm <- metadata

  # Extract relevant values
  meta_glm$y <- web@tableY[, which_dictionary[1]]
  meta_glm <- cbind(meta_glm, x = argonaut::getFeature(web@tableX.sft, which_dictionary[2]))

  all_terms <- attr(terms.formula(formula), "term.labels")
  all_subtypes <- colnames(meta_glm)[grep(x = colnames(meta_glm), pattern = "x\\.")]

  # All subtypes in one model
  vec_combined <- c(0, 1)
  # One model per subtype
  vec_full <- matrix(rep(c(0, 1), length(all_subtypes)), ncol = length(all_subtypes), byrow = F)
  # Differential association per subtype
  vec_out <- matrix(rep(c(0, 1), (2 * length(all_subtypes)) * length(all_terms)), ncol = length(all_subtypes), byrow = F)

  # paste together all subtypes with plusses between them:
  collapsed_subtypes <- paste(all_subtypes, collapse = " + ")

  internal_formula <- reformulate(
    termlabels = paste0("(", collapsed_subtypes, ") * 1 * ."),
    response = "y"
  )

  formula_full <- update.formula(formula, internal_formula)

  for (i in 1:length(all_subtypes)) {
    # paste together all subtypes with plusses between them:
    collapsed_subtypes_i <- paste(all_subtypes[i], collapse = " + ")

    internal_formula_i <- reformulate(
      termlabels = paste0("(", collapsed_subtypes_i, ") * 1 * ."),
      response = "y"
    )

    formula_i <- update.formula(formula, internal_formula_i)


    # fit linear model
    fit <- lm(formula = formula_i, data = meta_glm)

    # Calculate p-value for entire model
    fstat <- summary(fit)$fstatistic
    p <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

    vec_full[1, i] <- summary(fit)$r.squared
    vec_full[2, i] <- p
  }

  # fit linear model
  fit <- lm(formula = formula_full, data = meta_glm)

  # Calculate p-value for entire model
  fstat <- summary(fit)$fstatistic
  p <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

  vec_combined[1] <- summary(fit)$r.squared
  vec_combined[2] <- p

  # run ANOVA to determine impact of group on SLOPE of association.
  disj_fit <- anova(fit)

  # Calculate r.squared
  # get the residual sum of squares
  resid.ss <- sum(residuals(fit)^2)
  target_disj_interactions <- grep(paste(paste0("^", all_subtypes, ":"), collapse = "|"), row.names(disj_fit))
  disj.rsquared <- disj_fit[target_disj_interactions, 2] / (disj_fit[target_disj_interactions, 2] + resid.ss)

  vec_out[cbind(rep(x = c(1:length(all_terms) * 2) - 1, each = length(all_subtypes)), 1:(length(all_subtypes)))] <- disj.rsquared
  vec_out[cbind(rep(x = c(1:length(all_terms) * 2), each = length(all_subtypes)), 1:(length(all_subtypes)))] <- disj_fit[target_disj_interactions, 5]

  # For each subtype, run ANOVA to determine impact of group on STRENGTH of association.

  for (i in 1:length(all_subtypes)) {
    f_get_resid <- reformulate(
      termlabels = all_subtypes[i],
      response   = "y"
    )

    meta_glm$abs_resid <- abs(residuals(lm(formula = f_get_resid, data = meta_glm, na.action = na.exclude)))

    formula_emerg <- update.formula(
      new = reformulate(
        termlabels = paste0("(", paste(all_subtypes[-i], collapse = " + "), ") * 1 * ."),
        response   = "abs_resid"
      ), old = formula
    )

    emerg_fit <- lm(formula = formula_emerg, data = meta_glm)
    emerg_anova <- anova(emerg_fit)

    # Calculate r.squared
    resid_emerg.ss <- sum(residuals(emerg_fit)^2)

    target_emerg_interactions <- grep(paste(paste0("^", all_terms, "$"), collapse = "|"), row.names(emerg_anova))
    emerg.rsquared <- emerg_anova[target_emerg_interactions, 2] / (emerg_anova[target_emerg_interactions, 2] + resid_emerg.ss)

    vec_out[cbind(rep(x = 2 * length(all_terms) + c(1:length(all_terms) * 2) - 1, each = length(all_terms)), i)] <- emerg.rsquared
    vec_out[cbind(rep(x = 2 * length(all_terms) + c(1:length(all_terms) * 2), each = length(all_terms)), i)] <- emerg_anova[target_emerg_interactions, 5]
  }

  return(list(vec_combined = vec_combined, vec_full = vec_full, vec_out = vec_out))
}


#' Run differential correlation analysis for all interacting metabolites and functions in case of stratified output.
#' @description Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param stats_out A vector containing all statistical output from \code{model_picker}.
#' @param all_terms A vector containing all RHS terms of the supplied model formula, including interactions.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#'
collate_model_output_argonaut <- function(web, stats_out, all_terms) {
  # Create result container matrices. We take advantage of the fact that true interactions are coded as TRUE, which corresponds to 1,
  # automatically setting all non-canonical interactions as p = 1 and estimate = 0.
  dictionary <- web@dictionary
  strat_dict <- web@strat_dict
  colnames(strat_dict) <- paste(colnames(strat_dict), unlist(argonaut::apply_by(web@tableX.sft, 3, names)[1, ], use.names = F), sep = ".")

  out_rvals_tot <- dictionary
  out_pvals_tot <- !dictionary
  out_rvals <- strat_dict
  out_pvals <- !strat_dict
  out_disjrvals <- strat_dict
  out_disjpvals <- !strat_dict
  out_emergrvals <- strat_dict
  out_emergpvals <- !strat_dict

  stats_tot <- argonaut_parse_tot(stats_out = stats_out)

  out_rvals_tot[dictionary] <- stats_tot[, 1]
  out_pvals_tot[dictionary] <- stats_tot[, 2]

  # #Expand to stratified dimensions
  # out_rvals_tot = as.matrix(out_rvals_tot)[,rep(1:ncol(as.matrix(out_rvals_tot)), argonaut::apply_by(web@tableX.sft, 3, length)[1,])]
  # colnames(out_rvals_tot) <- colnames(strat_dict)
  # out_pvals_tot = as.matrix(out_pvals_tot)[,rep(1:ncol(as.matrix(out_pvals_tot)), argonaut::apply_by(web@tableX.sft, 3, length)[1,])]
  # colnames(out_pvals_tot) <- colnames(strat_dict)

  # Adjust for multiple comparisons
  out_qvals_tot <- out_pvals_tot
  out_qvals_tot[dictionary] <- NA

  out_tot_tale <- list(full = new("anansiTale",
    subject    = "model_full",
    type       = "r.squared",
    estimates  = out_rvals_tot,
    p.values   = out_pvals_tot,
    q.values   = out_qvals_tot
  ))

  stats_full <- argonaut_parse_full(stats_out = stats_out)
  stats_diff <- argonaut_parse_diff(stats_out = stats_out)

  out_rvals[strat_dict] <- stats_full[1, ]
  out_pvals[strat_dict] <- stats_full[2, ]

  # Adjust for multiple comparisons
  out_qvals <- out_pvals
  out_qvals[strat_dict] <- NA

  out_tale <- list(full = new("anansiTale",
    subject    = "model_full_by_subtype",
    type       = "r.squared",
    estimates  = out_rvals,
    p.values   = out_pvals,
    q.values   = out_qvals
  ))

  out_disjointed <- vector(mode = "list", length = length(all_terms))

  for (i in 1:length(all_terms)) {
    out_disjrvals[strat_dict] <- stats_diff[(i * 2) - 1, ]
    out_disjpvals[strat_dict] <- stats_diff[(i * 2), ]

    # Adjust for multiple comparisons
    out_disjqvals <- out_disjpvals
    out_disjqvals[strat_dict] <- NA


    out_disjointed[[i]] <- new("anansiTale",
      subject    = paste("model_disjointed", all_terms[i], sep = "_"),
      type       = "r.squared",
      estimates  = out_disjrvals,
      p.values   = out_disjpvals,
      q.values   = out_disjqvals
    )
  }
  names(out_disjointed) <- all_terms

  out_emergent <- vector(mode = "list", length = length(all_terms))

  for (i in 1:length(all_terms)) {
    out_emergrvals[strat_dict] <- stats_diff[(i * 2) + (i * 2) - 1, ]
    out_emergpvals[strat_dict] <- stats_diff[(i * 2) + (i * 2), ]

    # Adjust for multiple comparisons
    out_emergqvals <- out_emergpvals
    out_emergqvals[strat_dict] <- NA


    out_emergent[[i]] <- new("anansiTale",
      subject    = paste("model_emergent", all_terms[i], sep = "_"),
      type       = "r.squared",
      estimates  = out_emergrvals,
      p.values   = out_emergpvals,
      q.values   = out_emergqvals
    )
  }
  names(out_emergent) <- all_terms


  # Collect into nested list and return results
  return(list(
    modelfit = out_tot_tale,
    modelfit_sub = out_tale,
    disjointed = out_disjointed,
    emergent = out_emergent
  ))
}

#' Parse argonansi output
#' @description Typically, the main \code{anansi()} function will run this for you.
#' @param stats_out A vector containing all statistical output from \code{model_picker}.
#' @return a matrix, containing the parameters for the full combined models.
#'
argonaut_parse_tot <- function(stats_out) {
  do.call(rbind, lapply(stats_out, FUN = function(x) {
    x$vec_combined
  }))
}

#' Parse argonansi output
#' @description Typically, the main \code{anansi()} function will run this for you.
#' @param stats_out A vector containing all statistical output from \code{model_picker}.
#' @return a matrix, containing the parameters for the full individual models.
#'
argonaut_parse_full <- function(stats_out) {
  do.call(cbind, lapply(stats_out, FUN = function(x) {
    x$vec_full
  }))
}

#' Parse argonansi output
#' @description Typically, the main \code{anansi()} function will run this for you.
#' @param stats_out A vector containing all statistical output from \code{model_picker}.
#' @return a matrix, containing the parameters for the differential association models.
#'
argonaut_parse_diff <- function(stats_out) {
  do.call(cbind, lapply(stats_out, FUN = function(x) {
    x$vec_out
  }))
}

#' Collapse tableX for Argonaut data
#' @description Should not be run on its own. to be applied by \code{anansiDiffCor()}. Typically, the main \code{anansi()} function will handle this for you.
#' @param tableX.sft A stratifiedFeatureTable object.
#' @return a matrix with rows as samples and columns as ordered features per subtype.
#'
get_strat_X <- function(tableX.sft) {
  slice_list <- lapply(1:dim(tableX.sft)[2],
    FUN = function(f) {
      slic <- tableX.sft[, f, ][, !is.na(colSums(tableX.sft[, f, ]))]
      colnames(slic) <- paste(dimnames(tableX.sft)[[2]][f],
        colnames(slic),
        sep = "."
      )
      return(slic)
    }
  )
  slice_mat <- do.call(cbind, slice_list)
  slice_mat <- slice_mat[, order(colnames(slice_mat))]

  return(slice_mat)
}

#' Expand a matrix of dictionary dimensions to stratified dictionary dimensions.
#' @description Should not be run on its own.
#' @param x An anansiTale object with result matrices of dictionary dimensions.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @return An anansiTale object with result matrices with rows as tableX features and columns as ordered tableY features per subtype.
#'
argonaut_rep_to_strat <- function(x, web) {
  slot(x, "estimates") <- as.matrix(slot(x, "estimates"))[, rep(1:ncol(as.matrix(web@dictionary)), argonaut::apply_by(web@tableX.sft, 3, length)[1, ])]
  colnames(slot(x, "estimates")) <- paste(colnames(slot(x, "estimates")), unlist(argonaut::apply_by(web@tableX.sft, 3, names)[1, ], use.names = F), sep = ".")

  slot(x, "p.values") <- as.matrix(slot(x, "p.values"))[, rep(1:ncol(as.matrix(web@dictionary)), argonaut::apply_by(web@tableX.sft, 3, length)[1, ])]
  colnames(slot(x, "p.values")) <- paste(colnames(slot(x, "p.values")), unlist(argonaut::apply_by(web@tableX.sft, 3, names)[1, ], use.names = F), sep = ".")

  slot(x, "q.values") <- as.matrix(slot(x, "q.values"))[, rep(1:ncol(as.matrix(web@dictionary)), argonaut::apply_by(web@tableX.sft, 3, length)[1, ])]
  colnames(slot(x, "q.values")) <- paste(colnames(slot(x, "q.values")), unlist(argonaut::apply_by(web@tableX.sft, 3, names)[1, ], use.names = F), sep = ".")

  return(x)
}
#
#
# library(argonaut)
# library(tidyverse)
# set.seed(123)
#
# y = rnorm(100)
#
# x = board_argo(nsubtypes = 5,
#                nfeatures = 10,
#                nsamples = 100, p_missing = 0.3) %>%
#   as.sft()
# z = sample(LETTERS[1:3], 100, replace = T)
# u = sample(LETTERS[1:3], 100, replace = T)
#
# x = getFeature(x, 1)
# colnames(x) = paste0("x.", colnames(x))
# metadata = data.frame(y, z, u)
# f = ~ z *u
# metadata = cbind(metadata, x)
#
# formula = f
# meta_glm = metadata
#
# argonaut::apply_by(x, 3, names)
#
#
# temp = all_terms
#
# all_terms = c(temp)
# do.call(length, strsplit(all_terms, paste0(all_terms[!grepl(all_terms, pattern = ":")], collapse = "|"), ""))
# unlist(lapply(strsplit(all_terms, paste0(all_terms[!grepl(all_terms, pattern = "\\:")], collapse = "|"), ""), length))
