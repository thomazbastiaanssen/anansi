#' Calculate an association network
#' @description This is the main workspider function in the anansi package. It manages the individual functionalities of anansi, including correlation analysis, correlation by group and differential correlation.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor()} are "spearman" and "kendall".
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param adjust.method Method to adjust p-values for multiple comparisons. \code{adjust.method = "BH"} is the default value. See \code{p.adjust()} in the base R \code{stats} package.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @param diff_cor A boolean. Toggles whether to compute differential correlations. Default is \code{TRUE}.
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
#' KOs.exp = clr_lite(KOs)
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
#' results    = spinToWide(anansi_output = anansi_out)
#' }
#'
anansi = function(web, method = "pearson", groups = NULL, adjust.method = "BH", verbose = T, diff_cor = T){

  #generate anansiYarn output object
  outYarn = new("anansiYarn", input = new("anansiInput", web = web, groups = groups))

  assess = assessGroups(web = web, groups = groups, diff_cor = diff_cor)
  groups          = assess$groups
  diff_cor        = assess$diff_cor

  if(verbose){print("Running annotation-based correlations")}
  #initialize cor_output list object
  output = new("anansiOutput")

  output@cor_results        = anansiCorTestByGroup(web     = web,
                                                    method = method,
                                                    groups = groups,
                                                    adjust.method = adjust.method,
                                                    verbose = verbose)

  if(diff_cor){
  if(verbose){print("Fitting models for differential correlation testing")}
    output@model_results = anansiDiffCor(web = web, groups = groups, adjust.method = adjust.method)
  }
  outYarn@output = output

  return(outYarn)
}


#' Investigate validity of the groups argument
#' @description checks whether \code{groups} is missing, the correct length and suitable for correlations per groups and differential correlations.
#' If something looks off, \code{assessGroups} will do its best to salvage it and let you know something's up.
#' This is a helper function called by \code{anansi}.
#' @param web web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param diff_cor A boolean. Toggles whether to compute differential correlations. Default is \code{TRUE}.
#' @return a list including a modified \code{groups} and \code{diff_cor} argument.
#'
assessGroups <- function(web, groups, diff_cor = diff_cor){
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
