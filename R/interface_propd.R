#' Run propr for all interacting features within one sample.
#' Wraps around the propr function.
#' @description If the \code{groups} argument is suitable, will also run proportionality analysis per group. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param groups A categorical or continuous value necessary for differential proportionality Typically a state or treatment score.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total dataset and per group if applicable.
#' @importFrom propr propr
#' @seealso \code{\link{anansi}}
#'
wrap_propr <- function(web, groups, verbose = T) {
  # Determine all groups
  all_groups <- unique(groups)
  # If there are numbers here, we cannot do propr by group, so we'll substitute a groups called "All" for this part.
  if (!inherits(all_groups, "character")) {
    all_groups <- "All"
    groups <- rep("All", nrow(web@tableY))
  }

  # If that's all fine, we can determine the size of the output
  else if (length(all_groups) > 1) {
    all_groups <- c("All", all_groups)
  }

  # Generate container list of suitable length for all results
  out_list <- vector(mode = "list", length = length(all_groups))
  names(out_list) <- c(all_groups)

  # first run for all groups together
  out_list$All <- propr_by_group(web, groups = rep(T, nrow(web@tableY)), verbose = verbose)
  out_list$All@subject <- "All"


  # If verbose, verbalize.
  if (length(all_groups) > 1) {
    if (verbose) {
      print(paste(
        "Running proportionality for the following groups:",
        paste(all_groups, collapse = ", ")
      ))
    }

    # Assess proportions for subsets if applicable.
    # Skip 1 since that's taken by "All".
    for (i in 2:length(all_groups)) {
      out_by_group <- propr_by_group(web, groups = groups == all_groups[i], verbose = verbose)
      out_by_group@subject <- all_groups[i]

      out_list[[i]] <- out_by_group
    }
  }
  # Return results
  return(out_list)
}

#' Run propr for all interacting features such as metabolites or functions.
#' Wraps around the propr function.
#' @description If the \code{groups} argument is suitable, will also run proportionality analysis per group. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param groups A categorical or continuous value necessary for differential proportionality Typically a state or treatment score.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total dataset and per group if applicable.
#' @importFrom stats pt
#' @importFrom methods new
#' @importFrom propr propr getMatrix
#'
propr_by_group <- function(web, groups, verbose = T) {
  # Compute Rho
  pr <- propr::propr(web@tableY[groups, ], # rows as samples, like it should be
    metric = "rho", # or "phi", "phs", "cor", "vlr"
    ivar = NA
  ) # used by updateCutoffs

  # Extract Rho statistic
  r <- propr::getMatrix(pr) * web@dictionary

  # Compute t-statistics based on the n and the correlation coefficient
  n <- web@dictionary
  n[T] <- nrow(web@tableY[groups, ])
  t <- (r * sqrt(n - 2)) / sqrt(1 - r^2)

  # Compute p-values based on t and n.
  p <- 2 * (1 - pt(abs(t), (n - 2)))

  # Compute naive adjusted p-values
  q <- p
  q[web@dictionary] <- NA

  # Collate correlation coefficients, p-values and q-values into an anansiTale
  out <- new("anansiTale",
    subject    = "All",
    type       = "rho.values",
    estimates  = r,
    p.values   = p,
    q.values   = q
  )
  return(out)
}


# Disabling propd for now, leaving this to be sure:
#'
#' #' Run differential proportionality analysis using propr for all interacting features such as metabolites or functions.
#' #' Wraps around the propd function.
#' #' @description Can either take continuous or categorical data for \code{groups}. Typically, the main \code{anansi()} function will run this for you.
#' #' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' #' @param groups A categorical or continuous value necessary for differential proportionality Typically a state or treatment score.
#' #' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' #' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent proportionality and one for disjointed proportionality.
#' #' @importFrom propr propd
#' #' @importFrom future.apply future_apply
#' #'
#' differential_prop = function(web, groups, verbose = T){
#'
#'   #Call gather propd to compute thetas and p-values
#'   stats_out <- gather_propd(web = web, groups = groups, verbose = verbose)
#'
#'   #Create result container matrices. We take advantage of the fact that true interactions are coded as TRUE, which corresponds to 1,
#'   #automatically setting all non-canonical interactions as p = 1 and estimate = 0.
#'
#'   out_disjrvals  <- web@dictionary
#'   out_disjpvals  <- !web@dictionary
#'   out_emergrvals <- web@dictionary
#'   out_emergpvals <- !web@dictionary
#'
#'   out_disjrvals[web@dictionary]  <- stats_out$em_t[web@dictionary]
#'   out_disjpvals[web@dictionary]  <- stats_out$em_p[web@dictionary]
#'   out_emergrvals[web@dictionary] <- stats_out$dj_t[web@dictionary]
#'   out_emergpvals[web@dictionary] <- stats_out$dj_p[web@dictionary]
#'
#'
#'   out_disjqvals                  <- out_disjpvals
#'   out_disjqvals[web@dictionary]  <- NA
#'
#'   out_disjointed = new("anansiTale",
#'                        subject    = "model_disjointed",
#'                        type       = "theta_f",
#'                        estimates  = out_disjrvals,
#'                        p.values   = out_disjpvals,
#'                        q.values   = out_disjqvals)
#'
#'   out_emergqvals                 <- out_emergpvals
#'   out_emergqvals[web@dictionary] <- NA
#'
#'   out_emergent   = new("anansiTale",
#'                        subject    = "model_emergent",
#'                        type       = "theta_e",
#'                        estimates  = out_emergrvals,
#'                        p.values   = out_emergpvals,
#'                        q.values   = out_emergqvals)
#'
#'   #Collect into nested list and return results
#'   return(list(disjointed = out_disjointed,
#'               emergent   = out_emergent))
#' }
#'
#'
#' #' Run differential proportionality analysis using propr for all interacting features such as metabolites or functions.
#' #' Wraps around the propd function.
#' #' @description Can either take continuous or categorical data for \code{groups}. Typically, the main \code{anansi()} function will run this for you.
#' #' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' #' @param groups A categorical or continuous value necessary for differential proportionality Typically a state or treatment score.
#' #' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' #' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent proportionality and one for disjointed proportionality.
#' #' @importFrom propr propd
#' #' @importFrom future.apply future_apply
#' #'
#' wrap_propd = function(web, groups, verbose = T){
#'
#'   #Call gather propd to compute thetas and p-values
#'   stats_out <- gather_propd(web = web, groups = groups, verbose = verbose)
#'
#'   #Create result container matrices. We take advantage of the fact that true interactions are coded as TRUE, which corresponds to 1,
#'   #automatically setting all non-canonical interactions as p = 1 and estimate = 0.
#'
#'   out_disjrvals  <- web@dictionary
#'   out_disjpvals  <- !web@dictionary
#'   out_emergrvals <- web@dictionary
#'   out_emergpvals <- !web@dictionary
#'
#'   out_disjrvals[web@dictionary]  <- stats_out$em_t[web@dictionary]
#'   out_disjpvals[web@dictionary]  <- stats_out$em_p[web@dictionary]
#'   out_emergrvals[web@dictionary] <- stats_out$dj_t[web@dictionary]
#'   out_emergpvals[web@dictionary] <- stats_out$dj_p[web@dictionary]
#'
#'
#'   out_disjqvals                  <- out_disjpvals
#'   out_disjqvals[web@dictionary]  <- NA
#'
#'   out_disjointed = new("anansiTale",
#'                        subject    = "model_disjointed",
#'                        type       = "theta_f",
#'                        estimates  = out_disjrvals,
#'                        p.values   = out_disjpvals,
#'                        q.values   = out_disjqvals)
#'
#'   out_emergqvals                 <- out_emergpvals
#'   out_emergqvals[web@dictionary] <- NA
#'
#'   out_emergent   = new("anansiTale",
#'                        subject    = "model_emergent",
#'                        type       = "theta_e",
#'                        estimates  = out_emergrvals,
#'                        p.values   = out_emergpvals,
#'                        q.values   = out_emergqvals)
#'
#'   #Collect into nested list and return results
#'   return(list(disjointed = out_disjointed,
#'               emergent   = out_emergent))
#' }
#'
#' #' Run differential proportionality analysis using propr for all interacting features such as metabolites or functions.
#' #' Wraps around the propd function.
#' #' @description Can either take continuous or categorical data for \code{groups}. Typically, the main \code{anansi()} function will run this for you.
#' #' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' #' @param groups A categorical or continuous value necessary for differential proportionality Typically a state or treatment score.
#' #' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' #' @return a list with four matrices containing thetas and p-values for emergent and disjointed proportionality testing respectively.
#' #' @importFrom propr propd getMatrix setActive
#' #' @importFrom stats pf
#' #'
#' gather_propd = function(web, groups, verbose){
#'   # softmax to undo CLR.
#'   if(verbose){print("Applying softmax to CLR-transformed data...")}
#'   ct.list <- apply(web@tableY, 1, softmax, simplify = F)
#'   cts     <- do.call(rbind, ct.list)
#'
#'   #prepare for f-test
#'   K = length(unique(groups))
#'   N = length(groups)
#'
#'   pd <- propr::propd(counts = cts,
#'               group  = groups, # a vector of 2 or more groups
#'               alpha = NA, # whether to handle zeros
#'               weighted = FALSE # whether to weigh log-ratios
#'   )
#'
#'   pe <- propr::setActive(pd, what = "theta_e")
#'   pd <- propr::setActive(pd, what = "theta_f")
#'
#'   emerg_theta <- propr::getMatrix(pe)
#'   disj_theta  <- propr::getMatrix(pd)
#'   #emerg_theta = 1 - disj_theta
#'
#'   emerg_F <- theta_to_F(emerg_theta, N = N, K = K)
#'   disj_F  <- theta_to_F(disj_theta,  N = N, K = K)
#'
#'   emerg_p = pf(emerg_F, df1 = (K - 1), df2 = (N - K), lower.tail = FALSE)
#'   disj_p  = pf(disj_F,  df1 = (K - 1), df2 = (N - K), lower.tail = FALSE)
#'
#'
#'   return(list(em_t = emerg_theta,
#'               em_p = emerg_p,
#'               dj_t = disj_theta,
#'               dj_p = disj_p))
#' }
#'
#' #'Calculate F stat from theta. We will consider theta to be R^2-like as it is defined as the proportion of total variance that can be explained by the groups.
#' #'@param theta a theta statistic from propd
#' #'@param N number of observations
#' #'@param K number of groups
#' #'
#' theta_to_F = function(theta, N, K){
#' # F = (R^2 / (k-1)/
#' # (1 - R^2) / (n - k - 1)
#' #We will consider theta to be like R^2, so:
#'
#'
#'   return(   (( theta / (K -1) )/ ((1 - theta) / (N - K))))
#' }


# a = rnorm(36)
# b = sample(letters[1:3], 36, replace = T)
# anova(lm(a ~ b))
