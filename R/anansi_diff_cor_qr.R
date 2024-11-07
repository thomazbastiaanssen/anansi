#' Run differential correlation analysis for all interacting metabolites and functions.
#' @description Can either take continuous or categorical data for \code{groups}. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param formula A formula object. Used to assess differential associations.
#' @param reff A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.
#' @param modeltype A string, either "lm" or "lmer" depending on the type of model that should be ran.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals model.matrix.default terms.formula
#' @importFrom future.apply future_apply
#' @importFrom methods is
#'
anansiDiffCor = function(web, metadata, formula, reff, modeltype, verbose = T){
  #Create a matrix with row and column coordinates to cycle through the relevant comparisons in tableY and tableX.
  which_dictionary <- which(get_dict(web), arr.ind = T, useNames = F)
  lm.metadata      <- cbind(x = 1, metadata)

  all_terms  <- labels(terms.formula(formula))

  #make saturated model
  sat_model = update.formula(old = formula, ~ x * 1 * (.))
  if(verbose){cat(paste0("Fitting least-squares for following model:\n", paste0(as.character(sat_model), " ", collapse = "")))}

  #compute shape of model.matrix and initialize qr.mm
  contr <- make_contrasts(lm.metadata)
  base.mm <- suppressWarnings(
    model.matrix.default(sat_model, lm.metadata, contrasts.arg = contr)
  )

  #identify the columns of the variables that change based on X-variable
  x.fct  <- `dimnames<-`(attr(terms.formula(sat_model), "factors"), NULL)

  #TODO This could be important for argonaut support; consider modifying row id 1.
  x.assign   <- as.integer(which(x.fct[1,] == 1))
  all.assign <- attr(base.mm, "assign")

  x.vars     <- all.assign %in% x.assign
  #prevent the first round from removing anything.
  x.int      <- x.assign[-1]

  #set up F-tests
  #number of samples in full model
  n <- nrow(web@tableY)
  #number of parameters in full model
  sat_params <- colnames(base.mm)
  p_1 <- length(sat_params)

  #strip mm of attributes
  mm = matrix(base.mm, ncol = p_1)

  #Null models for general association of y and x: Remove all components with x.
  Y.TSS <- SS(fast.qr.resid(y = web@tableY,
                            x = mm[,!x.vars]))

  #prepare output
  df_mat <- dfmat(x.assign, x.int, all.assign, x.fct, n)

  modelfit  <- list(full = new("anansiTale",
                               subject    = "model_full",
                               type       = "r.squared",
                               df         = df_mat[,1],
                               estimates  = get_dict.double(web) * Y.TSS, #start with RSS0
                               F.values   = get_dict.double(web),
                               p.values   = !get_dict(web),
                               q.values   = !get_dict(web)))

  disjointed <- `names<-`(lapply(1:length(all_terms),
                                 function(x)
                                   new("anansiTale",
                                       subject    = paste("model_disjointed", all_terms[x], sep = "_"),
                                       type       = "r.squared",
                                       df         = df_mat[,x + 1],
                                       estimates  = get_dict.double(web), #start with RSS0
                                       F.values   = get_dict.double(web),
                                       p.values   = !get_dict(web),
                                       q.values   = !get_dict(web))), all_terms)

  emergent <- `names<-`(lapply(1:length(all_terms),
                                 function(x)
                                   new("anansiTale",
                                       subject    = paste("model_emergent", all_terms[x], sep = "_"),
                                       type       = "r.squared",
                                       df         = df_mat[,x + 1] + c(0, -1, -1/df_mat[1, x + 1]),
                                       estimates  = get_dict.double(web), #start with RSS0
                                       F.values   = get_dict.double(web),
                                       p.values   = !get_dict(web),
                                       q.values   = !get_dict(web))), all_terms)

#Compute R^2 for full model
modelfit$full@estimates <- 1 - (
  sapply(seq_len(NCOL(web@tableX)),
         function(x) R_full(y = x, mm, web, x.fct, x.vars)) /
    modelfit$full@estimates)

#compute all disjointed R^2 values, return to matrix with rows as values and columns as terms
for(t in seq_along(x.int)){

  i.disj <- index.disj(x = x.int[t], all.assign, x.fct);

  for(y in seq_len(NCOL(web@tableX))){
    #adjust the input model.matrix by multiplying the relevant columns by x
    qr.mm  <- mm
    qr.mm[,x.vars] <- mm[,x.vars] * web@tableX[,y]
    y.ind  <- get_dict(web)[,y]
    y.vals <- `dimnames<-`(web@tableY[,y.ind], NULL)

    #straight to web!!
    disjointed[[t]]@estimates[get_dict(web)[,y],y]  <- R_disj(y.vals, qr.mm, i.disj)

    emergent  [[t]]@estimates[get_dict(web)[,y],y]  <- R_emerg(y.vals, qr.mm, i.disj)
  }

}
#Add F and P statistics
modelfit   <- lapply(modelfit,   get_PF, d = get_dict(web))
disjointed <- lapply(disjointed, get_PF, d = get_dict(web))
emergent   <- lapply(emergent,   get_PF, d = get_dict(web))


return(list(
  modelfit = modelfit,
  disjointed = disjointed,
  emergent = emergent))
}

#' Get total sum of squares for residual
#' @noRd
#'
colwiseTSS <- function(x) apply(x, 2, function(x) sum((x - mean(x))^2))

#' fast resids
#' @noRd
#'
fast.qr.resid <- function(x,y) {
  y - crossprod(t(x), qr.coef(qr(x), y))
}

#' fast resids with sorted rolling
#' @noRd
#'
fast.qr.resid <- function(x,y) {
  y - crossprod(t(x), qr.coef(qr(x), y))
}


#' Get df1 & 2 for f-ratio
#' @noRd
#'
dfmat <- function(x.assign, x.int, all.assign, x.fct, n){
  df0   <- colSums(!sapply(x.assign, function(x) index.self.high(x, all.assign, x.fct)))
  df1   <- c(sum(index.self.high(x.assign[1], all.assign, x.fct)),
             sapply(x.int, function(x) sum(all.assign %in% x)))
  df2 = n - (df0 + df1)
  dfr = df2/df1
  return(rbind(df1,df2, dfr))
}

#' @noRd
#'
make_contrasts <- function(lm.metadata){
  contr.in <- NULL
f.names <- names(which(sapply(lm.metadata, function(x) is.character(x) || is.factor(x) || is.logical(x))))
if(length(f.names) > 0){
  contr.in <- `names<-`(rep("contr.sum", times = length(f.names)), f.names)
  contr.in[names(which(sapply(lm.metadata, is.ordered)))] <- "contr.poly"
}
return(as.list(contr.in))
}

#' Calculate sums of squares by column
#' @description Sums of Squared residuals, short for `colSums(x^2)`.
#' @noRd
#'
SS <- function(x) {
  dn <- dim(x)
  .colSums(x^2,dn[1L], dn[2L], FALSE)
}


#' Compute F and P statistic for \code{anansiTale} object.
#' @description Populate \code{anansiTale} object with F statistics.
#' @param object An \code{anansiTale} object.
#' @param d A binary adjacency matrix, corresponding to the relevant dictionary
#'
get_PF <- function(object, d){

  object@F.values[d] <- object@estimates[d] * object@df[3]
  object@p.values[d] <- pf(object@F.values[d],
                           df1 = object@df[1], df2 = object@df[2],
                           lower.tail = FALSE)

  return(object)
}

#' #' @noRd
#' #'
#' R_disj_emerg <- function(y, x, mm, web, x.vars, i.disj){
#'   #adjust the input model.matrix by multiplying the relevant columns by x
#'   qr.mm  <- mm
#'   qr.mm[,x.vars] <- mm[,x.vars] * web@tableX[,y]
#'   y.ind  <- get_dict(web)[,y]
#'   y.vals <- `dimnames<-`(web@tableY[,y.ind], NULL)
#'
#'   Rsq_d <- R_disj(y.vals, qr.mm, i.disj)
#'
#'   Rsq_e <- R_emerg(y.vals, qr.mm, i.disj)
#'
#'   return(cbind(Rsq_d,
#'                Rsq_e))
#'
#' }

#' @noRd
#'
  R_disj <- function(y.vals, qr.mm, i.disj){
# Disjointed
  mm.0 <- qr.mm[,i.disj[,1]]
  mm.1 <- qr.mm[,i.disj[,2]]


    #Cycle through dropping interactions with x, including higher order interactions.
    RSS_i0 <- SS(fast.qr.resid(y = y.vals,x = mm.0))

    #Cycle through returning interactions with x, but not higher order interactions.
    RSS_i1 <- SS(fast.qr.resid(y = y.vals, x = mm.1))

    return(1 - (RSS_i1 / RSS_i0))

  }

  #' @noRd
  #'
  R_emerg <- function(y.vals, qr.mm, i.disj){
    #Emergent
    e.ord   <- order(rowSums(qr.mm[,i.disj[,3], drop = FALSE]))
    e.mm    <- cbind(1, diff(qr.mm[e.ord,-1]))

    e.y     <- diff(y.vals[e.ord, drop = FALSE])

    #use diff to look at var
    e.mm.0  <- e.mm[,i.disj[,1]]
    e.mm.1  <- e.mm[,i.disj[,2]]

    #Cycle through dropping interactions with x, including higher order interactions.
    RSS_i0 <- SS(fast.qr.resid(y = e.y, x = e.mm.0))

    #Cycle through returning interactions with x, but not higher order interactions.
    RSS_i1 <- SS(fast.qr.resid(y = e.y, x = e.mm.1))

    return(1 - (RSS_i1 / RSS_i0))

  }


  #' @noRd
  #'
R_full <- function(y, mm, web, x.fct, x.vars){
  #adjust the input model.matrix by multiplying the relevant columns by x
  qr.mm  <- mm
  qr.mm[,x.vars] <- mm[,x.vars] * web@tableX[,y]
  y.ind  <- get_dict(web)[,y]
  y.vals <- `dimnames<-`(web@tableY[,y.ind], NULL)

  #Return H1 SSR
  SS(fast.qr.resid(y = y.vals, x = qr.mm))

}

#' @noRd
#'
index.self.high <- function(x, all.assign, x.fct){
  all.assign %in% which(apply(
    x.fct * x.fct[,x] == x.fct[,x],
    MARGIN =  2, all))
}

#' @noRd
#'
index.high <- function(x, all.assign, x.fct){

  all.assign %in% which(
    `[[<-`(apply(
    x.fct * x.fct[,x] == x.fct[,x],
    MARGIN =  2, all), subscript = x, FALSE)
    )
}

#' @noRd
#'
index.disj <- function(x, all.assign, x.fct){
  #first is full null,
  #second is parameter to investigate returned
  i0 <- !all.assign %in% which(apply(
    x.fct * x.fct[,x] == x.fct[,x],
    MARGIN =  2, all))
  ix <- all.assign %in% x
  cbind(i0, i1 = i0 | ix, ix)
}

