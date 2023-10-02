#' #' An S4 class to contain all metabolomics and functional input data as well as a dictionary to link them.
#' #' Can only be used with the argonaut library.
#' #'
#' #' @slot tableY A matrix of metabolomics data. Rows are samples and columns are features.
#' #' @slot tableX A stratified feature table containing features of interest. dim1 should be samples, dim2 should be features and dim3 should be subtypes.
#' #' @slot dictionary A binary adjacency matrix. Typically generated using the \code{weaveWebFromTables()} function.
#' #' @description argonansiWeb is the main container that will hold your input data thoughout the \code{anansi} pipeline.
#' #'
#' setClass("argonansiWeb",
#'          slots = c(
#'            tableY       = "matrix",
#'            tableX       = "matrix",
#'            tableX.sft   = "stratifiedFeatureTable",
#'            dictionary   = "matrix"
#'          )
#' )
#'
#'
#'
#' #' Create an argonansiWeb object from two 'omics tables and a dictionary
#' #' @description This function will take two tables of 'omics data, for example metabolomics and stratified functional microbiome data. It will also take a dictionary list, which is provided in this package.
#' #' @param tableY A table containing features of interest. Rows should be samples and columns should be features. The Y and X refer to the position of the features in a formula: Y ~ X.
#' #' @param tableX A stratified feature table containing features of interest. dim1 should be samples, dim2 should be features and dim3 should be subtypes.
#' #' @param dictionary A list that has feature names from tableY as names. Default is the dictionary provided in this package.
#' #' @param mode A character vector. Can be "interaction" or "membership". Toggles whether to link two datasets based on their interactions or based on shared group membership.
#' #' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' #' @param prune A boolean. Toggles whether to prune particularly large groups.
#' #' @param max_sds A numeric. Only relevant for prune == TRUE. How many SDs larger than the median a group can be before it is pruned.
#' #' For general use, we recommend sticking to that one. You can access the dictionary like this: \code{data(dictionary)}
#' #' @return an \code{anansiWeb} object. Web is used as input for most of the main workflow of anansi.
#' #' @export
#' weaveWebFromStratifiedTables = function(tableY, stratifiedTableX, dictionary = anansi::anansi_dic, verbose = T, mode = "interaction", prune = F, max_sds = 3){
#'
#'   stopifnot("stratifiedTableX needs to be a stratifiedFeatureTable object from the argonaut package." = "stratifiedFeatureTable" %in% class(stratifiedTableX))
#'
#'   tableX <- argonaut::apply_by(X = stratifiedTableX, MARGIN = 3, mean)
#'
#'   if(length(dictionary) == 1){
#'     if(dictionary == "none"){
#'       if(verbose){print("No dictionary provided, preparing for all vs all analysis. ")}
#'       dictionary <- mock_dictionary(tableY = tableY, tableX = tableX)
#'     }
#'   }
#'
#'   stopifnot("the mode argument needs to be interaction or membership." = mode %in% c("interaction", "membership"))
#'   #For conventional use, table Y should be metabolites and table X functions.
#'
#'   # The two 'table' matrices MUST have row
#'   # and column names that are unique, and
#'   # look like the following:
#'   #
#'   #               f1  f2  f3  f4  f5  f6
#'   #   sample1      0   0   2   0   0   1
#'   #   sample2     20   8  12   5  19  26
#'   #   sample3      3   0   2   0   0   0
#'   #       ... many more rows ...
#'   #
#'
#'   if(prune){
#'     dictionary = prune_dictionary_for_exclusivity(dict_list = dictionary,
#'                                                   max_sds = max_sds, verbose = verbose)
#'   }
#'
#'   #create binary adjacency matrix first
#'   dictionary = makeAdjacencyMatrix(tableY = tableY, dict_list = dictionary,
#'                                    verbose = verbose, mode = mode)
#'
#'   #If we're looking at single data set, don't do associations with yourself. Set diagonal to FALSE.
#'   if(identical(tableX, tableY)){
#'     diag(dictionary) <- F
#'   }
#'
#'   #Check if input tables have the same names and the same length.
#'   if(!identical(row.names(tableY), row.names(tableX)))
#'   {warning("The row names of tableY and tableX do not correspond. Please make sure they are in the same order.")}
#'   if(nrow(tableY) != nrow(tableX))
#'   {stop("tableY and tableX do not have the same number of rows/observations.")}
#'
#'   available_tableY = sort(intersect(colnames(tableY), rownames(dictionary)))
#'   available_tableX = sort(intersect(colnames(tableX), colnames(dictionary)))
#'
#'   if(verbose){
#'     print(paste(length(available_tableY), "were matched between table 1 and the columns of the adjacency matrix"))
#'     print(paste(length(available_tableX), "were matched between table 2 and the rows of the adjacency matrix"))
#'   }
#'   #Select the relevant part of the adjacency matrix
#'   dictionary = dictionary[available_tableY, available_tableX]
#'
#'   #Ensure there are no features that never interact
#'   dictionary = dictionary[rowSums(dictionary) > 0,colSums(dictionary) > 0]
#'
#'   #use the dictionary to clean input tables
#'   tableY = tableY[,rownames(dictionary)]
#'   tableX = tableX[,colnames(dictionary)]
#'   stratifiedTableX <- stratifiedTableX[,colnames(dictionary),]
#'
#'   #Return an argonansiWeb object with four slots: typically metabolites, functions, stratified functions and adjacency matrix
#'   return(new("argonansiWeb",
#'              tableY     = as.matrix(tableY),
#'              tableX     = as.matrix(tableX),
#'              tableX.sft = stratifiedTableX,
#'              dictionary = as.matrix(dictionary)))
#' }
#'
#'
#'
#' #' Run differential correlation analysis for all interacting metabolites and functions.
#' #' @description Should not be run on its own. to be applied by \code{anansiDiffCor()}. Typically, the main \code{anansi()} function will handle this for you.
#' #' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' #' @param which_dictionary A matrix derived from calling \code{which(web@dictionary, arr.ind = T)}. It contains coordinates for the relevant measurements to be compared.
#' #' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' #' @param formula A formula object. Used to assess differential associations.
#' #' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' #' @importFrom stats anova lm pf residuals na.exclude
#' #'
#' glm_argonaut_calc_diff_cor <- function(web, which_dictionary, metadata, formula){
#'
#'   meta_glm = metadata
#'
#'   # Extract relevant values
#'   meta_glm$y = web@tableY[,which_dictionary[1]]
#'   meta_glm = cbind(meta_glm, x = argonaut::getFeature(web@tableX.sft[,which_dictionary[2],]))
#'
#'   all_terms    <- attr(terms(formula), "term.labels")
#'   all_subtypes <- colnames(meta_glm)[grep(x = colnames(meta_glm), pattern = "x\\.")]
#'   vec_out = rep(c(0, 1), 1 + 2 * (length(all_subtypes) + length(all_subtypes) * length(all_terms)))
#'
#'   #paste together all subtypes with plusses between them:
#'   collapsed_subtypes = paste(all_subtypes, collapse = " + ")
#'
#'   internal_formula = reformulate(
#'     termlabels = paste0("(",collapsed_subtypes ,") * 1 * ."),
#'     response    = "y")
#'
#'   formula_full = update.formula(formula, internal_formula)
#'
#'   # fit linear model
#'   fit      <- lm(formula = formula_full, data = meta_glm)
#'
#'   # Calculate p-value for entire model
#'   fstat    <- summary(fit)$fstatistic
#'   p        <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
#'
#'   vec_out[1] <- summary(fit)$r.squared
#'   vec_out[2] <- p
#'
#'   #run ANOVA to determine impact of group on SLOPE of association.
#'   disj_fit        <- anova(fit)
#'
#'   #Calculate r.squared
#'   #get the residual sum of squares
#'   resid.ss                 <- sum(residuals(fit)^2)
#'   target_disj_interactions <- grep("^x:",row.names(disj_fit))
#'   disj.rsquared            <- disj_fit[target_disj_interactions,2] / (disj_fit[target_disj_interactions,2] + resid.ss)
#'
#'   vec_out[1 + (1:(length(all_subtypes) + length(all_subtypes) * length(all_terms))*2)]  <- disj.rsquared
#'   vec_out[2 + (1:(length(all_subtypes) + length(all_subtypes) * length(all_terms))*2)]  <- disj_fit[target_disj_interactions,5]
#'
#'   #run ANOVA to determine impact of group on STRENGTH of association.
#'   formula_emerg = update.formula(abs_resid ~ ., formula)
#'
#'   meta_glm$abs_resid      <- abs(residuals(lm(y ~ x, data = meta_glm, na.action = na.exclude)))
#'
#'   emerg_fit      <- lm(formula = formula_emerg, data = meta_glm)
#'   emerg_anova    <- anova(emerg_fit)
#'
#'   #Calculate r.squared
#'   resid_emerg.ss <- sum(residuals(emerg_fit)^2)
#'   target_emerg_interactions <- 1:length(all_terms)
#'
#'   emerg.rsquared <- emerg_anova[target_emerg_interactions,2] / (emerg_anova[target_emerg_interactions,2] + resid_emerg.ss)
#'
#'   vec_out[1 + 2 * length(all_terms) + (1:length(all_terms))*2] <- emerg.rsquared
#'   vec_out[2 + 2 * length(all_terms) + (1:length(all_terms))*2] <- emerg_anova[target_emerg_interactions,5]
#'
#'   return(vec_out)
#' }
#'
#'
#'
#' library(argonaut)
#' library(tidyverse)
#' set.seed(123)
#'
#' y = rnorm(100)
#'
#' x = board_argo(nsubtypes = 5,
#'                nfeatures = 10,
#'                nsamples = 100, p_missing = 0.3) %>%
#'   as.sft()
#' z = sample(LETTERS[1:3], 100, replace = T)
#'
#' df.2 = data.frame(z = z, y = y, x = getFeature(x, "Philosophy"))
#' meta_glm = df.2
