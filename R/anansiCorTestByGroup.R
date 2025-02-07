#' Run correlations for all interacting metabolites and functions.
#' @description If the \code{groups} argument is suitable, will also run correlation analysis per group. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param group.vec A character vector denoting group membership. Typically a state or treatment score.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total dataset and per group if applicable.
#' @seealso \code{\link{anansi}}
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
    group.bool = rep(TRUE, NROW(get_tableY(web))), verbose = verbose
  )

  out_list$All@subject <- "All"

  if (!is.null(all_groups)) {
    # If verbose, verbalize.
    if (verbose) {
      message(paste(
        "Running correlations for the following groups:",
        paste(all_groups, collapse = ", ")
      ))
    }
    for (i in seq_len(length(all_groups))) {
      out_by_group <- anansiCorPvalue(web, group.bool = group.vec == all_groups[i], verbose = verbose)
      out_by_group@subject <- all_groups[i]

      out_list[[i + 1]] <- out_by_group
    }
  }

  # Return results
  return(out_list)
}

#' Compute r-statistics for each featureY-featureX pair in the dictionary.
#' Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param group.bool A categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return An \code{anansiTale} result object.
#' @seealso \code{\link{anansi}} \cr \code{\link{anansiCorTestByGroup}}
#' @importFrom stats pt
#' @importFrom methods new
#'
anansiCorPvalue <- function(web, group.bool, verbose = verbose) {
  # Compute correlation coefficients
  r <- anansiCor(web = web, group.bool = group.bool)

  # Compute t-statistics based on the n and the correlation coefficient
  n <- sum(group.bool)
  t <- (r * sqrt(n - 2)) / sqrt(1 - r^2)

  # Compute p-values based on t and n.
  p <- 2 * (1 - pt(abs(t), (n - 2)))

  # Compute naive adjusted p-values
  q <- p
  q[web@dictionary] <- NA

  # Collate correlation coefficients, p-values and q-values into an anansiTale
  out <- new("anansiTale",
    subject    = "All",
    type       = "r.values",
    estimates  = r,
    p.values   = p,
    q.values   = q
  )
  return(out)
}

#' Compute r-statistics for each featureY-featureX pair in the dictionary.
#' Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param group.bool A boolean vector used to select which samples should be included in the correlations.
#' @seealso \code{\link{anansi}} \cr \code{\link{anansiCorTestByGroup}}
#' @return A matrix of r-statistics.
#' @importFrom stats cor
#'
anansiCor <- function(web, group.bool) {
  # Run correlations on subsections of your data
  cors <- cor(
    x = get_tableY(web)[group.bool, ],
    y = get_tableX(web)[group.bool, ],
    method = "pearson", use = "pairwise.complete.obs"
  )
  # set non-canonical correlations to zero using the binary adjacency matrix.
  return(cors * get_dict(web))
}
