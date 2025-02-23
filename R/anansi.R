#' Calculate an association network
#' @description This is the main workspider function in the anansi package. It
#' manages the individual functionalities of anansi, including correlation
#' analysis, correlation by group and differential correlation.
#' @param web An \code{anansiWeb} object, containing two tables with 'omics data
#' and a dictionary that links them. See \code{weaveWebFromTables()} for how to
#' weave a web.
#' @param metadata A vector or data.frame of categorical or continuous value
#' necessary for differential correlations. Typically a state or treatment.
#' If no argument provided, anansi will let you know and still to regular
#' correlations according to your dictionary.
#' @param groups A vector of the column names of categorical values in the
#' metadata object to control which groups should be assessed for simple
#' correlations.
#' If no argument provided, anansi will let you know and still to regular
#' correlations according to your dictionary.
#' @param formula A formula object. Used to assess differential associations.
#' @param adjust.method Method to adjust p-values for multiple comparisons.
#' \code{adjust.method = "BH"} is the default value. See \code{p.adjust()} in
#' the base R \code{stats} package.
#' @param verbose A boolean. Toggles whether to print diagnostic information
#' while running. Useful for debugging errors on large datasets.
#' @param ignore_dictionary A boolean. Default is FALSE. If set to TRUE, regular
#' all vs all associations will be tested regardless of the dictionary.
#' @param return.format \code{Character scalar}. Should be one of \code{"table"}
#' , \code{"list"}, or \code{"raw"}. Should the output of \code{\link{anansi}}
#' respectively be a wide `data.frame` of results, a list containing the results
#' and input, or a list of raw output (used for testing purposes).
#' convenient use. (Default: \code{"table"})
#' @return A list of lists containing correlation coefficients, p-values and
#' q-values for all operations.
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
#'   groups = "Legend",
#'   metadata = FMT_metadata,
#'   adjust.method = "BH",
#'   verbose = TRUE
#' )
#'
#' # To recreate the long plot:
#' library(tidyr)
#' library(ggplot2)
#'
#'anansiLong <- anansi_out |>
#'    pivot_longer(starts_with("All") | contains("FMT")) |>
#'    separate_wider_delim(
#'        name, delim = "_", names = c("cor_group", "cor_param")
#'            ) |>
#'            pivot_wider(names_from = cor_param, values_from = value)
#'
#' # Now it's ready to be plugged into ggplot2, though let's filter a bit more.
#'
#' # Only consider interactions where the entire model fits well enough.
#' anansiLong <- anansiLong[anansiLong$full_q.values < 0.1, ]
#'
#'
#'
#' ggplot(
#'   data = anansiLong,
#'   aes(
#'     x = r.values,
#'     y = feature_X,
#'     fill = cor_group,
#'     alpha = disjointed_Legend_p.values < 0.05
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
                   adjust.method = "BH", verbose = TRUE,
                   ignore_dictionary = FALSE, return.format = "table") {

  # generate anansiYarn input object
  input <- prepInput(
    web = web, formula = formula, groups = groups,
    metadata = metadata, verbose = verbose
  )
  int.terms <- input@int.terms; groups <- input@groups; n.grps <- input@n.grps;
  group.id <- input@group.id; errorterm <- input@error.term;
  sat_model <- input@lm.formula

  out.list <- vector(
    "list", length = 1 + n.grps + (2 * length(int.terms))
    )

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

  out.list[seq_len(n.grps)] <- call_groupwise(
    web, groups,
    metadata, verbose
  )
  # Sort out metadata formatting for differential association testing
  meta.frame <- model.frame(formula = sat_model, cbind(x = 1, metadata))

  out.list[n.grps +
           seq_len(1 + (2 * length(int.terms)))] <- anansiDiffCor(
             web, sat_model, errorterm, int.terms, meta.frame, verbose)

  if(return.format != "raw") {
    results <- result.df(out.list, get_dict(web))
    results <- anansi.p.adjust(results, adjust.method)
    attr(results, "group_terms") <- named_group_list(group.id, groups, metadata)
    attr(results, "model_terms") <- named_term_list(int.terms, metadata)
  }

  switch(return.format,
         "table" = return(results),
         "list"  = return(list(results, input = input)),
         "raw"   = return(out.list))
}


#' Manages group-wise association calls
#' @description If the \code{groups} argument is suitable, will also run correlation analysis per group. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @noRd
#'
call_groupwise <- function(web, groups, metadata, verbose) {
  if(is.null(groups)) {group.vec <- NULL} else {
  group.vec <- apply(metadata[,groups, drop = FALSE], 1, paste, collapse = "_")
  }

  if (is(web, "argonansiWeb")) {
    return(anansiCorTestByGroup(
      web = new("anansiWeb",
        tableY     = as.matrix(get_tableY(web)),
        tableX     = as.matrix(get_tableX(web)),
        dictionary = as.matrix(web@strat_dict)
      ), group.vec = group.vec, verbose = verbose
    ))
  }

  return(
    anansiCorTestByGroup(web, group.vec, verbose)
    )
}


#' Assess formula, trim metadata and prepare output for anansi workflow
#' @description Initialize `anansiInput` component of output. Should not be
#' called by user.
#' @noRd
#'
prepInput <- function(web, formula, groups, metadata, verbose) {
  raw_terms <- terms.formula(formula, "Error", data = metadata)
  indErr <- attr(raw_terms, "specials")$Error

  groups <- check_groups(groups, raw_terms, indErr, metadata, verbose)

  all_terms <- if (is.null(indErr)) {
    labels(raw_terms)
  } else {
    labels(raw_terms)[-indErr]
  }

  sat_model <- make_saturated_model(formula, raw_terms, indErr, verbose)
  error.term <- if (is.null(indErr)) {
    NULL
  } else {
    deparse1(attr(raw_terms, "variables")[[1L + indErr]][[2L]], backtick = TRUE)
  }
    input = new("anansiInput",
      web = web,
      lm.formula = sat_model,
      error.term = error.term,
      int.terms = all_terms,
      groups = groups[[1]],
      n.grps = groups[[2]],
      group.id = c("All",unique(apply(metadata[,groups[[1]], drop = FALSE],
                                       1, paste, collapse = "_"))),
      metadata = `row.names<-.data.frame`(metadata, NULL)
    )

  return(input)
}

#' Check group argument.
#' @noRd
#' @importFrom stats as.formula update.formula
#'
check_groups <- function(groups, raw_terms, indErr, metadata, verbose) {
  # If user input, check it
  if (!is.null(groups)) {

    missing_groups <- !groups %in% colnames(metadata)
    stopifnot(
      "Grouping variable(s) not recognised. Please check input and labels. " =
        !any(missing_groups)
      )
    n.groups <- 1 + length(unique(do.call(paste0, c(metadata[groups]))))
    return(list(groups, n.groups))
  }

  # If no input, look for categorical variables
  ind_o1 <- attr(raw_terms, "order") == 1
  if (!is.null(indErr)) {
    ind_o1[indErr] <- FALSE
  }
  if (!any(ind_o1)) {
    if(verbose){message("No grouping variable found for groupwise correlations. ")}
    return(list(NULL, 1))
  }
  groups <- labels(raw_terms)[ind_o1]

  sub_meta <- metadata[, groups, drop = FALSE]
  sub_meta <-
    sub_meta[, unname(apply( sub_meta, 2, is.categorical)), drop = FALSE]

  if (NCOL(sub_meta) == 0) {
    if(verbose) message("No grouping variable found for groupwise correlations.")
    return(list(NULL, 1))
  }
  n.groups <- 1 + length(unique(do.call(paste0, c(sub_meta))))

  return(list(groups, n.groups))
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

#' Is character, factor or ordered factor
#' @description wrapper around
#' \code{is.character(x) || is.factor(x) || is.ordered(x)}
#' @param x an object to be evaluated as being categorical
#' @returns a boolean.
#'
is.categorical <- function(x) is.character(x) || is.factor(x) || is.ordered(x)

#' Return levels or 'numeric'
#' @noRd
#'
lvs_or_num <- function(x) {
  if(is.character(x)) return(unique(x))
  return("numeric")}

#' Provide name and levels for categories for group terms
#' @description convenience function to generate group terms output
#' @param g character vector of all 'All' followed by all unique terms
#' @param t character vector of term(s)
#' @param m data.frame of provided metadata
#' @returns a named list of character vectors, where names are terms and
#' containing vectors are unique levels or 'numeric' if not categorical. First
#' Entry is named 'All' and contains all levels.
#'
named_group_list <- function(g, t, m) c(list(All = g), lapply(m[t], lvs_or_num))

#' Provide name and levels for categories
#' @description convenience function to generate output
#' @param t character vector of term(s)
#' @param m data.frame of provided metadata
#' @returns a named list of character vectors, where names are terms and
#' containing vectors are unique levels or 'numeric' if not categorical.
#'
named_term_list <- function(t, m) {
  order_one <- t %in% colnames(m)
  f_order   <- t[ order_one]
  h_order   <- t[!order_one]

  f_list <- lapply(m[f_order], lvs_or_num)
  h_list <- NULL
  if(length(h_order) != 0) {
    m_num <- !unlist(lapply(m, is.categorical))
    m[m_num] <- "numeric"
    h_terms <- strsplit(h_order, split = ":", fixed = TRUE)
    h_list  <- lapply(h_terms,
                      function(x) unique(do.call(paste, c(m[x], sep = "_"))))
    names(h_list) <- h_order
  }

  return(c(f_list, h_list))
}
