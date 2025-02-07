#' Handle FDR methods for anansi. Can also be used on anansi output to recalculate FDR.
#' @param x an AnansiYarn object, the main output from anansi().
#' @param method The p-value adjustment method. See ?p.adjust.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @importFrom stats p.adjust
#' @importFrom methods slot slot<-
#' @noRd
#'
anansiAdjustP <- function(x, method = "BH", verbose = TRUE) {
  # first check if the settings make sense
  stopifnot("Input needs to be an anansiYarn object, the output of the anansi() function" = is(x, "anansiYarn"))

  # verbalise the plan if verbose
  if (verbose) {
    if (method %in% c("BH", "fdr")) {
      message("Adjusting p-values using Benjamini & Hochberg's procedure.")
    } else {
      message("Adjusting p-values...")
    }
  }

  # extract dictionary for convenience
  dictionary <- yarn.dic(x)
  if (is(yarn.web(x), "argonansiWeb")) {
    dictionary <- x@input@web@strat_dict
  }

  # Take inventory of where to find p-values.
  pval_df <- data.frame(
    model = c(
      rep("cor_results", length(x@output@cor_results)),
      rep("model_results", length(x@output@model_results))
    ),
    group = c(
      names(x@output@cor_results),
      names(x@output@model_results)
    )
  )

  if (is(x@input@web, "argonansiWeb")) {
    pval_df <- pval_df[pval_df$group != "modelfit.full", ]
    p <- slot(slot(x@output, "model_results")[["modelfit.full"]], "p.values")

      # adjust p-values directly into the relevant slot
      slot(slot(x@output, "model_results")[["modelfit.full"]], "q.values")[x@input@web@dictionary] <-
        p.adjust(p[x@input@web@dictionary], method = method)
  }

  # for each of those sources of p-values, do:
  for (i in seq_len(NROW(pval_df))) {
    # Extract p-values
    p <- slot(slot(x@output, pval_df[i, 1])[[pval_df[i, 2]]], "p.values")

    # adjust p-values directly into the relevant slot
    slot(slot(x@output, pval_df[i, 1])[[pval_df[i, 2]]], "q.values")[dictionary] <-
      p.adjust(p[dictionary], method = method)

  }
  return(x)
}

