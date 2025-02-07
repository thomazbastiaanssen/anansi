#' Run differential correlation analysis for all interacting metabolites and functions.
#' @description Can either take continuous or categorical data for \code{groups}. Typically, the main \code{anansi()} function will run this for you.
#' @param yarn An \code{anansiYarn} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals model.matrix.default terms.formula
#' @importFrom future.apply future_apply
#' @importFrom methods is
#'
anansiDiffCor <- function(yarn, metadata, verbose = TRUE) {
  # Create a matrix with row and column coordinates to cycle through the
  # relevant comparisons in tableY and tableX.
  which_dictionary <- which(yarn.dic(yarn), arr.ind = TRUE, useNames = FALSE)

  # Compute shape of model.matrix and initialize qr.mm
  base.mm <- build.mm(sat_model = yarn.f(yarn), metadata)

  # identify the columns of the variables that change based on X-variable
  x.fct <- get_x.fct(yarn.f(yarn), errorterm = yarn.e(yarn))

  # TODO This could be important for argonaut support; consider modifying row id 1.
  x.assign <- as.integer(which(x.fct[1, ] == 1))
  all.assign <- attr(base.mm, "assign")
  x.vars <- all.assign %in% x.assign
  x.int <- x.assign[-1]

  # set up F-tests
  # number of samples in full model
  n <- nrow(yarn.tY(yarn))

  # strip mm of attributes
  mm <- matrix(base.mm, ncol = NCOL(base.mm))

  # Null models for general association of y and x: Remove all components with x.
  Y.TSS <- SS(fast.qr.resid(
    y = yarn.tY(yarn),
    x = mm[, !x.vars]
  ))

  # prepare output
  df_mat <- dfmat(x.assign, x.int, all.assign, x.fct, n)

  modelfit <- list(full = new("anansiTale",
    subject    = "model_full",
    type       = "r.squared",
    df         = df_mat[, 1],
    estimates  = yarn.dic.double(yarn) * Y.TSS, # start with RSS0
    F.values   = yarn.dic.double(yarn),
    p.values   = !yarn.dic(yarn),
    q.values   = !yarn.dic(yarn)
  ))

  disjointed <- `names<-`(lapply(
    seq_len(length(yarn.int(yarn))),
    function(x) {
      new("anansiTale",
        subject    = paste("model_disjointed", yarn.int(yarn)[x], sep = "_"),
        type       = "r.squared",
        df         = df_mat[, x + 1],
        estimates  = yarn.dic.double(yarn), # start with RSS0
        F.values   = yarn.dic.double(yarn),
        p.values   = !yarn.dic(yarn),
        q.values   = !yarn.dic(yarn)
      )
    }
  ), yarn.int(yarn))

  emergent <- `names<-`(lapply(
    seq_len(length(yarn.int(yarn))),
    function(x) {
      new("anansiTale",
        subject    = paste("model_emergent", yarn.int(yarn)[x], sep = "_"),
        type       = "r.squared",
        df         = df_mat[, x + 1] + c(0, -1, -1 / df_mat[1, x + 1]),
        estimates  = yarn.dic.double(yarn), # start with RSS0
        F.values   = yarn.dic.double(yarn),
        p.values   = !yarn.dic(yarn),
        q.values   = !yarn.dic(yarn)
      )
    }
  ), yarn.int(yarn))

  # Compute R^2 for full model
  modelfit$full@estimates <- 1 - (
    sapply(
      seq_len(NCOL(yarn.tX(yarn))),
      function(x) R_full(y = x, mm, yarn, x.fct, x.vars)
    ) /
      modelfit$full@estimates)

  # compute all disjointed R^2 values, return to matrix with rows as values and columns as terms
  for (t in seq_along(x.int)) {
    i.disj <- index.disj(x = x.int[t], all.assign, x.fct)

    for (y in seq_len(NCOL(yarn.tX(yarn)))) {
      # adjust the input model.matrix by multiplying the relevant columns by x
      qr.mm <- mm
      qr.mm[, x.vars] <- mm[, x.vars] * yarn.tX(yarn)[, y]
      y.ind <- yarn.dic(yarn)[, y]
      y.vals <- `dimnames<-`(yarn.tY(yarn)[, y.ind], NULL)

      # straight to web!!
      disjointed[[t]]@estimates[yarn.dic(yarn)[, y], y] <- R_disj(y.vals, qr.mm, i.disj)

      emergent[[t]]@estimates[yarn.dic(yarn)[, y], y] <- R_emerg(y.vals, qr.mm, i.disj)
    }
  }
  # Add F and P statistics
  modelfit <- lapply(modelfit, get_PF, d = yarn.dic(yarn))
  disjointed <- lapply(disjointed, get_PF, d = yarn.dic(yarn))
  emergent <- lapply(emergent, get_PF, d = yarn.dic(yarn))


  return(list(
    modelfit = modelfit,
    disjointed = disjointed,
    emergent = emergent
  ))
}

#' Get total sum of squares for residual
#' @noRd
#'
colwiseTSS <- function(x) apply(x, 2, function(x) sum((x - mean(x))^2))

#' fast resids
#' @noRd
#'
fast.qr.resid <- function(x, y) {
  y - crossprod(t(x), qr.coef(qr(x), y))
}

#' fast resids with sorted rolling
#' @noRd
#'
fast.qr.resid <- function(x, y) {
  y - crossprod(t(x), qr.coef(qr(x), y))
}


#' Get df1 & 2 for f-ratio
#' @noRd
#'
dfmat <- function(x.assign, x.int, all.assign, x.fct, n) {
  df0 <- colSums(
    !vapply(x.assign,
      FUN.VALUE = vector("logical", length = length(all.assign)),
      function(x) index.self.high(x, all.assign, x.fct)
    )
  )
  df1 <- c(
    sum(index.self.high(x.assign[1], all.assign, x.fct)),
    vapply(x.int, FUN.VALUE = 2, function(x) sum(all.assign %in% x))
  )
  df2 <- n - (df0 + df1)
  dfr <- df2 / df1
  return(rbind(df1, df2, dfr))
}

#' @noRd
#'
make_contrasts <- function(metadata) {
  contr.in <- NULL
  f.names <- colnames(metadata)[unlist(lapply(
    metadata,
    function(x) (is.character(x) || is.factor(x) || is.logical(x))
  ))]
  if (length(f.names) > 0) {
    contr.in <- `names<-`(rep("contr.sum", times = length(f.names)), f.names)
    contr.in[names(which(vapply(metadata, FUN.VALUE = FALSE, is.ordered)))] <- "contr.poly"
  }
  return(as.list(contr.in))
}

#' @noRd
#'
build.mm <- function(sat_model, metadata) {
  contr <- make_contrasts(metadata)

  return(model.matrix.default(
    sat_model, metadata,
    contrasts.arg = contr
  ))
}


#' Calculate sums of squares by column
#' @description Sums of Squared residuals, short for `colSums(x^2)`.
#' @noRd
#'
SS <- function(x) {
  dn <- dim(x)
  .colSums(x^2, dn[1L], dn[2L], FALSE)
}

#' Transform to odds: `x / (1-x)`.
#' @noRd
#'
oddify <- function(x) x / (1 - x)

#' Compute F and P statistic for \code{anansiTale} object.
#' @description Populate \code{anansiTale} object with F statistics.
#' @param object An \code{anansiTale} object.
#' @param d A binary adjacency matrix, corresponding to the relevant dictionary
#' @noRd
#'
get_PF <- function(object, d) {
  object@F.values[d] <- oddify(object@estimates[d]) * object@df[3]
  object@p.values[d] <- pf(object@F.values[d],
    df1 = object@df[1], df2 = object@df[2],
    lower.tail = FALSE
  )

  return(object)
}


#' @noRd
#'
R_disj <- function(y.vals, qr.mm, i.disj) {
  # Disjointed
  mm.0 <- qr.mm[, i.disj[, 1]]
  mm.1 <- qr.mm[, i.disj[, 2]]


  # Cycle through dropping interactions with x, including higher order interactions.
  RSS_i0 <- SS(fast.qr.resid(y = y.vals, x = mm.0))

  # Cycle through returning interactions with x, but not higher order interactions.
  RSS_i1 <- SS(fast.qr.resid(y = y.vals, x = mm.1))

  return(1 - (RSS_i1 / RSS_i0))
}

#' @noRd
#'
R_emerg <- function(y.vals, qr.mm, i.disj) {
  # Emergent
  e.ord <- order(rowSums(qr.mm[, i.disj[, 3], drop = FALSE]))
  e.mm <- cbind(1, diff(qr.mm[e.ord, -1]))

  e.y <- diff(y.vals[e.ord, drop = FALSE])

  # use diff to look at var
  e.mm.0 <- e.mm[, i.disj[, 1]]
  e.mm.1 <- e.mm[, i.disj[, 2]]

  # Cycle through dropping interactions with x, including higher order interactions.
  RSS_i0 <- SS(fast.qr.resid(y = e.y, x = e.mm.0))

  # Cycle through returning interactions with x, but not higher order interactions.
  RSS_i1 <- SS(fast.qr.resid(y = e.y, x = e.mm.1))

  return(1 - (RSS_i1 / RSS_i0))
}


#' @noRd
#'
R_full <- function(y, mm, yarn, x.fct, x.vars) {
  # adjust the input model.matrix by multiplying the relevant columns by x
  qr.mm <- mm
  qr.mm[, x.vars] <- mm[, x.vars] * yarn.tX(yarn)[, y]
  y.ind <- yarn.dic(yarn)[, y]
  y.vals <- `dimnames<-`(yarn.tY(yarn)[, y.ind], NULL)

  # Return H1 SSR
  SS(fast.qr.resid(y = y.vals, x = qr.mm))
}

#' @noRd
#'
index.self.high <- function(x, all.assign, x.fct) {
  all.assign %in% which(apply(
    x.fct * x.fct[, x] == x.fct[, x],
    MARGIN = 2, all
  ))
}

#' @noRd
#'
index.high <- function(x, all.assign, x.fct) {
  all.assign %in% which(
    `[[<-`(apply(
      x.fct * x.fct[, x] == x.fct[, x],
      MARGIN = 2, all
    ), subscript = x, FALSE)
  )
}

#' @noRd
#'
index.disj <- function(x, all.assign, x.fct) {
  # first is full null,
  # second is parameter to investigate returned
  i0 <- !all.assign %in% which(apply(
    x.fct * x.fct[, x] == x.fct[, x],
    MARGIN = 2, all
  ))
  ix <- all.assign %in% x
  cbind(i0, i1 = i0 | ix, ix)
}

#' Generate x.fct table, deal with \code{Error} terms.
#' @noRd
#'
get_x.fct <- function(sat_model, errorterm) {
  x.fct <- attr(terms.formula(sat_model), "factors")

  if (!is.null(errorterm)) {
    x.fct <- x.fct[rownames(x.fct) != errorterm, ]
  }

  return(`dimnames<-`(x.fct, NULL))
}

#' Generate table, to compute contrasts, basis for base.mm
#' @noRd
#'
subset_metadata <- function(metadata, keep, raw_terms, indErr) {
  if (!is.null(indErr)) {
    all_terms <- c(
      all_terms, deparse1(attr(raw_terms, "variables")[[1L + indErr]][[2L]],
        backtick = TRUE
      )
    )
  }
  return(cbind(x = 1, metadata[, colnames(metadata) %in% all_terms, drop = FALSE]))
}
