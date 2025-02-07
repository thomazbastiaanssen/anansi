#' Calculate an association network
#' @description This is the main workspider function in the anansi package. It manages the individual functionalities of anansi, including correlation analysis, correlation by group and differential correlation.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param groups A vector of the column names of categorical values in the metadata object to control which groups should be assessed for simple correlations. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param formula A formula object. Used to assess differential associations.
#' @param adjust.method Method to adjust p-values for multiple comparisons. \code{adjust.method = "BH"} is the default value. See \code{p.adjust()} in the base R \code{stats} package.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @param ignore_dictionary A boolean. Default is FALSE. If set to TRUE, regular all vs all associations will be tested regardless of the dictionary.
#' @return A list of lists containing correlation coefficients, p-values and q-values for all operations.
#' @importFrom stats model.frame
#' @export
#' @examples
#' # Load example data:
#'
#' data(dictionary)
#' data(FMT_data)
#'
#' # Clean and prepare the example data.
#' # In the example dataset, the metabolites are already cleaned.
#'
#' KOs <- floor(FMT_KOs)
#' KOs <- apply(KOs, c(1, 2), function(x) as.numeric(as.character(x)))
#' KOs <- KOs[apply(KOs == 0, 1, sum) <= (ncol(KOs) * 0.90), ]
#'
#' KOs <- KOs[row.names(KOs) %in% sort(unique(unlist(anansi_dic))), ]
#'
#' # CLR-transform.
#'
#' KOs.exp <- clr_c(KOs)
#'
#' # Make sure that columns are features and rows are samples.
#'
#' t1 <- t(FMT_metab)
#' t2 <- t(KOs.exp)
#'
#' # Run anansi pipeline.
#'
#' web <- weaveWebFromTables(
#'   tableY = t1,
#'   tableX = t2,
#'   dictionary = anansi_dic
#' )
#'
#' anansi_out <- anansi(
#'   web = web,
#'   formula = ~Legend,
#'   groups = FMT_metadata$Legend,
#'   metadata = FMT_metadata,
#'   adjust.method = "BH",
#'   verbose = TRUE
#' )
#'
#' results <- spinToWide(
#'   anansi_output = anansi_out, translate = TRUE,
#'   Y_translation = anansi::cpd_translation,
#'   X_translation = anansi::KO_translation
#' )
#'
#' # To recreate the long plot:
#' library(ggplot2)
#'
#' anansiLong <- spinToLong(
#'   anansi_output = anansi_out, translate = TRUE,
#'   Y_translation = anansi::cpd_translation,
#'   X_translation = anansi::KO_translation
#' )
#'
#' # Now it's ready to be plugged into ggplot2, though let's clean up a bit more.
#'
#' # Only consider interactions where the entire model fits well enough.
#' anansiLong <- anansiLong[anansiLong$model_full_q.values < 0.1, ]
#'
#'
#'
#' ggplot(
#'   data = anansiLong,
#'   aes(
#'     x = r.values,
#'     y = feature_X,
#'     fill = type,
#'     alpha = model_disjointed_Legend_p.values < 0.05
#'   )
#' ) +
#'
#'   # Make a vertical dashed red line at x = 0
#'   geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
#'
#'   # Points show  raw correlation coefficients
#'   geom_point(shape = 21, size = 3) +
#'
#'   # facet per compound
#'   ggforce::facet_col(~feature_Y, space = "free", scales = "free_y") +
#'
#'   # fix the scales, labels, theme and other layout
#'   scale_y_discrete(limits = rev, position = "right") +
#'   scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1 / 3)) +
#'   scale_fill_manual(values = c(
#'     "Young yFMT" = "#2166ac",
#'     "Aged oFMT" = "#b2182b",
#'     "Aged yFMT" = "#ef8a62",
#'     "All" = "gray"
#'   )) +
#'   theme_bw() +
#'   ylab("") +
#'   xlab("Pearson's rho")
#'
#' # See also ?spinToPlots
#'
anansi <- function(web, formula = ~1, groups = NULL, metadata,
                   adjust.method = "BH", verbose = TRUE, ignore_dictionary = FALSE) {

  # generate anansiYarn output object
  outYarn <- prepYarn(web = web, formula = formula, groups = groups,
                      metadata = metadata, verbose = verbose)
  # sort out metadata
  metadata <- model.frame(formula = yarn.f(outYarn), cbind(x = 1, metadata))

  if (ignore_dictionary) {
    if (verbose) {
      message("Dictionary will be ignored. Running all vs all associations.")
    }
    # set dictionary to all TRUE
    web@dictionary <- web@dictionary == web@dictionary
    if (is(web, "argonansiWeb")) {
      web@strat_dict <- web@strat_dict == web@strat_dict
    }
  }

  # initialize cor_output list object
  output <- new("anansiOutput")

  output@cor_results <- call_groupwise(
    web = web, groups = yarn.grp(outYarn),
    verbose = verbose
  )

  output@model_results <- unlist(anansiDiffCor(
      yarn = outYarn, metadata = metadata,
      verbose = verbose
    ))
  outYarn@output <- output

  # FDR
  outYarn <- anansiAdjustP(x = outYarn,
                           method = adjust.method, verbose = verbose)

  return(outYarn)
}


#' Manages group-wise association calls
#' @description If the \code{groups} argument is suitable, will also run correlation analysis per group. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @noRd
#'
call_groupwise <- function(web, groups, verbose) {

    if (is(web, "argonansiWeb")) {
      return(anansiCorTestByGroup(
        web = new("anansiWeb",
          tableY     = as.matrix(get_tableY(web)),
          tableX     = as.matrix(get_tableX(web)),
          dictionary = as.matrix(web@strat_dict)
        ), groups = groups, verbose = verbose
      ))
    }

    return(anansiCorTestByGroup(web = web, groups = groups, verbose = verbose))
  }


#' Assess formula, trim metadata and prepare output for anansi workflow
#' @description Initialize anansiYarn output Should not be called by user.
#' @noRd
#'
prepYarn <- function(web, formula, groups, metadata, verbose){

  raw_terms <- terms.formula(formula, "Error", data = metadata)
  indErr  <- attr(raw_terms, "specials")$Error

  groups <- check_groups(groups, raw_terms, indErr, metadata)

  all_terms <- if(is.null(indErr)) {labels(raw_terms)
    } else {labels(raw_terms)[-indErr]}

  sat_model <- make_saturated_model(formula, raw_terms, indErr, verbose)
  error.term <- if( is.null(indErr))  {NULL} else {
    deparse1(attr(raw_terms, "variables")[[1L + indErr]][[2L]], backtick = TRUE)}
  outYarn <- new("anansiYarn",
                 input = new("anansiInput",
                             web = web,
                             lm.formula = sat_model,
                             error.term = error.term,
                             int.terms  = all_terms,
                             groups = groups)
                 )
  return(outYarn)
}

#' Check group argument.
#' @noRd
#' @importFrom stats as.formula update.formula
#'
check_groups <- function(groups, raw_terms, indErr, metadata){
  if(!is.null(groups)) return(unname(groups))

  order_1 <- attr(raw_terms, "order") == 1
  if( !is.null(indErr) ) {order_1[indErr] <- FALSE}
  if( length(order_1) == !is.null(indErr) ) {return(NULL)}

  sub_meta <- metadata[, labels(raw_terms)[order_1], drop = FALSE]
  sub_meta <- sub_meta[,unname(apply(
    sub_meta, 2, function(x) {
      is.character(x) || is.factor(x) || is.ordered(x) })), drop = FALSE]

  if(NCOL(sub_meta) == 0) {return(groups = NULL)}
  groups <- apply(sub_meta, 1, paste, collapse = "_")
  return(unname(groups))
}

#' Prepare saturated model, deal with \code{Error} terms.
#' @noRd
#' @importFrom stats as.formula update.formula
#'
make_saturated_model <- function(formula, raw_terms, indErr, verbose) {
  # Simple case; No random intercept; make regular saturated model
  if (is.null(indErr)) {
    sat_model <- update.formula(old = formula, ~ x * 1 * (.))
    if (verbose) {
      message(
        paste0(
          "Fitting least-squares for following model:\n",
          paste0(as.character(sat_model), " ", collapse = "")
        )
      )
    }
    return(sat_model)
  }

  # Case with repeated measures:
 stopifnot("Only one Error() term allowed; more detected." = length(indErr) < 2)
  errorterm <- attr(raw_terms, "variables")[[1L + indErr]]
  sat_model <- update.formula(old = formula, new = as.formula(
    paste(
      "~",
      deparse1(errorterm[[2L]], backtick = TRUE),
      "+ x * 1 * (. -",
      deparse1(errorterm, backtick = TRUE),
      ")"
    ),
    env = environment(formula)
  ))
  if (verbose) {
    message(paste0(
      "Fitting least-squares for following model:\n",
      paste0(as.character(update.formula(old = formula, new = as.formula(
        paste("~ x * 1 * (. -", deparse1(errorterm, backtick = TRUE), ")"),
        env = environment(formula)
      ))), " ", collapse = ""),
      "\nwith '",
      deparse1(errorterm[[2L]], backtick = TRUE),
      "' as random intercept."
    ))
  }
  return(sat_model)
}
