#' Calculate an association network
#' @description This is the main workspider function in the anansi package. It manages the individual functionalities of anansi, including correlation analysis, correlation by group and differential correlation.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor()} are "spearman" and "kendall".
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param group A vector of the column names of categorical values in the metadata object to control which groups should be assessed for simple correlations. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param formula A formula object. Used to assess differential associations.
#' @param reff A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.
#' @param modeltype A string, either "lm" or "lmer" depending on the type of model that should be ran, or "propr" in the case of within-composition associations..
#' @param adjust.method Method to adjust p-values for multiple comparisons. \code{adjust.method = "BH"} is the default value. See \code{p.adjust()} in the base R \code{stats} package.
#' @param resampling A boolean. For p-value adjustment. Toggles the resampling of p-values to help reduce the influence of dependence of p-values. Will take more time on large datasets.
#' @param locality A boolean. For p-value adjustment. Toggles whether to prefer sampling from p-values from a comparison that shares an x or y feature. In a nutshell, considers the local p-value landscape when more important when correcting for FDR.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @param diff_cor A boolean. Toggles whether to compute differential correlations. Default is \code{TRUE}.
#' @param ignore_dictionary A boolean. Default is FALSE. If set to TRUE, regular all vs all associations will be tested regardless of the dictionary.
#' @return A list of lists containing correlation coefficients, p-values and q-values for all operations.
#' @export
#' @examples
#' \dontrun{
#' #Load example data:
#'
#' data(dictionary)
#' data(FMT_data)
#'
#' #Clean and prepare the example data.
#' #In the example dataset, the metabolites are already cleaned.
#'
#' KOs   <- floor(FMT_KOs)
#' KOs   <- apply(KOs,c(1,2),function(x) as.numeric(as.character(x)))
#' KOs   <- KOs[apply(KOs == 0, 1, sum) <= (ncol(KOs) * 0.90), ]
#'
#' KOs   <- KOs[row.names(KOs) %in% sort(unique(unlist(anansi_dic))),]
#'
#' #CLR-transform.
#'
#' KOs.exp = clr_c(KOs)
#'
#' #Make sure that columns are features and rows are samples.
#'
#' t1 = t(FMT_metab)
#' t2 = t(KOs.exp)
#'
#' #Run anansi pipeline.
#'
#' web        = weaveWebFromTables(tableY     = t1,
#'                                 tableX     = t2,
#'                                 dictionary = anansi_dic)
#'
#' anansi_out = anansi(web     = web,
#'                     method  = "pearson",
#'                     groups  = FMT_metadata$Legend,
#'                     adjust.method = "BH",
#'                     verbose = TRUE)
#'
#' results    = spinToWide(anansi_output = anansi_out, translate = T,
#'                         Y_translation = anansi::cpd_translation,
#'                         X_translation = anansi::KO_translation)
#'
#' #To recreate the long plot:
#' library(ggplot2)
#'
#' anansiLong <- spinToLong(anansi_output = anansi_out, translate = T,
#'                          Y_translation = anansi::cpd_translation,
#'                          X_translation = anansi::KO_translation)
#'
#' #Now it's ready to be plugged into ggplot2, though let's clean up a bit more.
#'
#' #Only consider interactions where the entire model fits well enough.
#' anansiLong <- anansiLong[anansiLong$model_full_q.values < 0.1,]
#'
#'
#'
#' ggplot(data = anansiLong,
#'        aes(x      = r.values,
#'            y      = feature_X,
#'            fill   = type,
#'            alpha  = model_disjointed_p.values < 0.05)) +
#'
#' #Make a vertical dashed red line at x = 0
#' geom_vline(xintercept = 0, linetype = "dashed", colour = "red")+
#'
#' #Points show  raw correlation coefficients
#' geom_point(shape = 21, size = 3) +
#'
#' #facet per compound
#' ggforce::facet_col(~feature_Y, space = "free", scales = "free_y") +
#'
#' #fix the scales, labels, theme and other layout
#' scale_y_discrete(limits = rev, position = "right") +
#' scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1/3)) +
#' scale_fill_manual(values = c("Young yFMT" = "#2166ac",
#'                              "Aged oFMT"  = "#b2182b",
#'                              "Aged yFMT"  = "#ef8a62",
#'                              "All"        = "gray"))+
#' theme_bw() +
#' ylab("") +
#' xlab("Pearson's rho")
#'
#' #See also ?spinToPlots
#' }
#'
anansi = function(web, method = "pearson", groups = NULL,  metadata = NULL, formula = ~1,
                  adjust.method = "BH", modeltype = "lm", resampling = F, locality = F,
                  reff = NULL, verbose = T, diff_cor = T, ignore_dictionary = F){

  if(is.vector(metadata)){
    metadata = data.frame(metadata = metadata)
  }

  #If there is a formula that is not empty:
  if(!identical(all.vars(formula), character(0))){
    unique_factors = all.vars(update.formula(formula, 0~.))

    if(is.null(groups)){
      groups = unique_factors
    }

    metadata_factors = metadata[,groups]
    groups = do.call(paste, c(data.frame(metadata_factors), sep = "_"))

    if(length(c(groups)) == 0){
      groups <- NULL
    }

  }



  if(ignore_dictionary){
    if(verbose){print("Dictionary will be ignored. Running all vs all associations.")}
    #set dictionary to all TRUE
    web@dictionary <- web@dictionary == web@dictionary
  }

  #assess validity of input
  assess_res      = assessAnansiCall(web = web, groups = groups, diff_cor = diff_cor,
                                     modeltype = modeltype, reff = reff, method = method)

  #adjust model call if appropriate
  groups          = assess_res$groups
  diff_cor        = assess_res$diff_cor
  reff            = assess_res$reff
  method          = assess_res$method

  #generate anansiYarn output object
  outYarn = new("anansiYarn", input = new("anansiInput", web = web, groups = groups, reff = reff))

  if(verbose){print("Running annotation-based correlations")}

  #initialize cor_output list object
  output = new("anansiOutput")

  output@cor_results     = call_groupwise(web     = web, method  = method,
                                          groups  = groups, verbose = verbose)

  if(diff_cor){
  if(verbose){print("Fitting models for differential correlation testing")
    print(paste("Model type:", modeltype, sep = ""))}
    output@model_results = anansiDiffCor(web = web, groups = groups, reff = reff,
                                         modeltype = modeltype, verbose = verbose)
  }
  outYarn@output = output

  #FDR
  outYarn <- anansiAdjustP(x = outYarn, method = adjust.method, resampling = resampling, locality = locality, verbose = verbose)

  return(outYarn)
}


#' Investigate validity of the Anansi Call
#' @description Calls both assessGroups and assessModelType.
#' This is a helper function called by \code{anansi}.
#' @param web web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param groups A vector of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param diff_cor A boolean. Toggles whether to compute differential correlations. Default is \code{TRUE}.
#' @param reff A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.
#' @param modeltype A string, either "lm" or "lmer" depending on the type of model that should be ran.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor()} are "spearman" and "kendall".
#' @return a list including a modified \code{groups} and \code{diff_cor} argument.
#'
assessAnansiCall <- function(web, groups, diff_cor = diff_cor, modeltype, reff, method){
  #create output list
  assess_out <- list(groups = groups, diff_cor = diff_cor, reff = reff, method = method)

  #assess validity of input
  assessG         = assessGroups(web = web, groups = groups, diff_cor = diff_cor)

  assess_out$groups          = assessG$groups
  assess_out$diff_cor        = assessG$diff_cor

  assessM         = assessModelType(modeltype = modeltype, reff = reff, diff_cor = diff_cor)

  assess_out$diff_cor        = assessM$diff_cor
  assess_out$reff            = assessM$reff

  if(identical(web@tableY, web@tableX)){
    assess_out$method = "propr"
  }

  return(assess_out)
}

#' Investigate validity of the groups argument
#' @description checks whether \code{groups} is missing, the correct length and suitable for correlations per groups and differential correlations.
#' If something looks off, \code{assessGroups} will do its best to salvage it and let you know something's up.
#' This is a helper function called by \code{anansi}.
#' @param web web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param diff_cor A boolean. Toggles whether to compute differential correlations. Default is \code{TRUE}.
#' @return a list including a modified \code{groups} and \code{diff_cor} argument.
#'
assessGroups <- function(web, groups, diff_cor = diff_cor){

  #first let's catch any trouble with factors.
if (inherits(groups, "factor")){
  groups <- as.character(groups)
  }
  #Three main options for data coming in:
  #1) groups is empty/NA
  #2) groups is numerical for diff_cor only
  #3) groups is categorical for diff_cor and cor_by_group

#Check if groups has the same length as the number of observations.
if(!is.null(groups) & length(groups) != nrow(web@tableY)){
  warning("The `groups` argument has a different length than your input tables. \nSomething is likely wrong, please check your input. ")
  diff_cor = F
  groups = rep("All", nrow(web@tableY))
}

#check if groups is missing or broken
if(diff_cor){

  if(is.null(groups) | any(is.na(groups))){
    message("Please be aware that the `groups` argument is missing or contains NAs. \nAnansi will proceed without differential association testing")
    diff_cor = F
    groups = rep("All", nrow(web@tableY))
  }
}

#In the case that groups is categorical, check if there are enough observations per group.
if(inherits(groups, "character")){
  if(!all(table(groups) > 3)) {
    warning("The `groups` argument is categorical, but not all groups have at least three observations.
              Correlations per group cannot be done. Please check your groups. ")
    diff_cor = F
    groups = rep("All", nrow(web@tableY))
  }

}

return(list(diff_cor = diff_cor,
            groups   = groups))
}

#' Investigate validity of the model call
#' @description checks the model call is legitimate
#' @param reff A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.
#' @param modeltype A string, either "lm" or "lmer" depending on the type of model that should be ran.
#' @param diff_cor A boolean. Toggles whether to compute differential correlations. Default is \code{TRUE}.
#'
assessModelType <- function(modeltype, reff, diff_cor){
  if(diff_cor){
    if(!modeltype %in% c("lmer", "lm")){
      warning("modeltype was not recognised. Needs to be exaclty `lm`, or `lmer`. Disabling differential association analysis. ")
      diff_cor = FALSE
    }
    if(modeltype == "lmer" & is.null(reff)){
      warning("lmer needs a categorical value for reff depicting group membership. Disabling differential association analysis. ")
      diff_cor = FALSE
    }
    if(modeltype != "lmer" & !is.null(reff)){
      warning("Only lmer can take a random effect. Either activate lmer or remove the reff argument. Disabling differential association analysis. ")
      diff_cor = FALSE
    }
  }
  if(is.null(reff)){reff = NA}
  return(list(diff_cor = diff_cor,
              reff     = reff))
}

#' Manages group-wise association calls
#' @description If the \code{groups} argument is suitable, will also run correlation analysis per group. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param method Correlation method.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#'
call_groupwise <- function(web, method, groups, verbose){
  if(method %in% c("pearson", "kendall", "spearman")){
    return(anansiCorTestByGroup(web = web, method = method, groups = groups, verbose = verbose))
  }

  if(method == "propr"){
    stopifnot("Two different tables detected, propr is meant for within-composition association testing." = identical(web@tableY, web@tableX))
    return(wrap_propr(web = web, groups = groups, verbose = verbose))
  }

}
