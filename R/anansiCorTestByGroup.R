#' Manages group-wise association calls
#' @description If the `groups` argument is suitable, will also run correlation analysis per group. Typically, the main `anansi()` function will run this for you.
#' @param web An `AnansiWeb` object, containing two tables with omics data and a dictionary that links them. See `weaveWebFromTables()` for how to weave a web.
#' @param groups A categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @noRd
#'
call_groupwise <- function(web, groups, metadata, verbose) {
  if(is.null(groups)) {group.vec <- NULL} else {
    group.vec <- apply(metadata[,groups, drop = FALSE], 1, paste, collapse = "_")
  }
  return(
    anansiCorTestByGroup(web, group.vec, verbose)
  )
}

#' Run correlations for all interacting metabolites and functions.
#' @description If the `groups` argument is suitable, will also run correlation analysis per group. Typically, the main `anansi()` function will run this for you.
#' @param web An `AnansiWeb` object, containing two tables with omics data and a dictionary that links them. See `weaveWebFromTables()` for how to weave a web.
#' @param group.vec A character vector denoting group membership. Typically a state or treatment score.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of `anansiTale` result objects, one for the total dataset and per group if applicable.
#' @seealso [anansi()]
#'
anansiCorTestByGroup <- function(web, group.vec, verbose = TRUE) {
  # Determine all groups
  all_groups <- unique(group.vec)

  # Generate container list of suitable length for all results
  out_list <- vector(mode = "list", length = 1 + length(all_groups))
  names(out_list) <- c("All", all_groups)

  # first run for all groups together
  out_list$All <- anansiCorPvalue(
    web,
    group.bool = rep(TRUE, NROW(web@tableY)), verbose
  )

  if (!is.null(all_groups)) {
    # If verbose, verbalize.
    if (verbose) {
      message(paste(
        "Running correlations for the following groups:\n",
        paste(all_groups, collapse = ", ")
      ))
    }
    for (i in seq_along(all_groups)) {
      out_by_group <- anansiCorPvalue(
        web, group.bool = group.vec == all_groups[i], verbose
        )
      out_by_group@subject <- all_groups[i]
      out_list[[i + 1]] <- out_by_group
    }
  }

  # Return results
  return(out_list)
}

#' Compute r-statistics for each featureY-featureX pair in the dictionary.
#' Typically, the main `anansi()` function will run this for you.
#' @param web An `AnansiWeb` object, containing two tables with omics data and a dictionary that links them. See `weaveWebFromTables()` for how to weave a web.
#' @param group.bool A categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return An `anansiTale` result object.
#' @seealso [anansi()] \cr [anansiCorTestByGroup()]
#' @importFrom stats pt
#' @importFrom methods new
#'
anansiCorPvalue <- function(web, group.bool, verbose) {
  # Compute correlation coefficients
  r <- anansiCor(web = web, group.bool = group.bool)

  # Compute p-value through t-statistic, based on n and correlation coefficient.
  n <- sum(group.bool)
  t <- abs((r * sqrt(n - 2)) / sqrt(1 - r^2))
  p <- 2 * (1 - pt(t, (n - 2)))

  # Collate correlation coefficients, p-values and q-values into an anansiTale
  out <- new("anansiTale",
    subject     = "All",
    type        = "r.values",
    estimates   = r,
    df          = n - 2,
    t.values    = t,
    p.values    = p
  )
  return(out)
}

#' Compute r-statistics for each featureY-featureX pair in the dictionary.
#' Typically, the main `anansi()` function will run this for you.
#' @param web An `AnansiWeb` object, containing two tables with omics data and a dictionary that links them. See `weaveWebFromTables()` for how to weave a web.
#' @param group.bool A boolean vector used to select which samples should be included in the correlations.
#' @seealso [anansi()] \cr [anansiCorTestByGroup()]
#' @return A matrix of r-statistics.
#' @importFrom stats cor
#'
anansiCor <- function(web, group.bool) {
  # Run correlations on subsections of your data
  cors <- cor(
    x = web@tableY[group.bool, ],
    y = web@tableX[group.bool, ],
    method = "pearson", use = "pairwise.complete.obs"
  )
  cors[! Matrix::as.matrix(web@dictionary) ] <- 0
  # set non-canonical correlations to zero using the binary adjacency matrix.
  return(cors)
}
