#' Run differential correlation analysis for all interacting metabolites and functions.
#' @description Can either take continuous or categorical data for \code{groups}. Typically, the main \code{anansi()} function will run this for you.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param formula A formula object. Used to assess differential associations.
#' @param reff A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.
#' @param modeltype A string, either "lm" or "lmer" depending on the type of model that should be ran.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return a list of \code{anansiTale} result objects, one for the total model, one for emergent correlations and one for disjointed correlations.
#' @importFrom stats anova lm pf residuals terms
#' @importFrom future.apply future_apply
#'
anansiDiffCor2 = function(web, metadata, formula, reff, modeltype, verbose = T){
  #Create a matrix with row and column coordinates to cycle through the relevant comparisons in tableY and tableX.
  which_dictionary <- which(web@dictionary, arr.ind = T, useNames = F)
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
  x.fct  <- `dimnames<-`(attr(terms(sat_model), "factors"), NULL)

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
                               df         = df_mat[,1]),
                               estimates  = web@dictionary * Y.TSS, #start with RSS0
                               F.values   = web@dictionary,
                               p.values   = !web@dictionary,
                               q.values   = !web@dictionary))

  disjointed <- `names<-`(lapply(1:length(all_terms),
                                 function(x)
                                   new("anansiTale",
                                       subject    = paste("model_disjointed", all_terms[x], sep = "_"),
                                       type       = "r.squared",
                                       df         = df_mat[,x + 1],
                                       estimates  = web@dictionary, #start with RSS0
                                       F.values   = web@dictionary,
                                       p.values   = !web@dictionary,
                                       q.values   = !web@dictionary)), all_terms)

  emergent <- `names<-`(lapply(1:length(all_terms),
                                 function(x)
                                   new("anansiTale",
                                       subject    = paste("model_emergent", all_terms[x], sep = "_"),
                                       type       = "r.squared",
                                       df         = df_mat[,x + 1] + c(0, -1, -1/df_mat[1, x + 1]),
                                       estimates  = web@dictionary, #start with RSS0
                                       F.values   = web@dictionary,
                                       p.values   = !web@dictionary,
                                       q.values   = !web@dictionary)), all_terms)

#Compute R^2 for full model
modelfit$full@estimates <- 1 - (
  sapply(seq_len(NCOL(web@tableX)),
         function(x) R_full(y = x, mm, web, x.fct, x.vars)) /
    modelfit$full@estimates)

#compute all disjointed R^2 values, return to matrix with rows a s values and columns as terms
for(t in seq_along(x.int)){

  i.disj <- index.disj(x = x.int[t], all.assign, x.fct);

  for(y in seq_len(NCOL(web@tableX))){
    #adjust the input model.matrix by multiplying the relevant columns by x
    qr.mm  <- mm
    qr.mm[,x.vars] <- mm[,x.vars] * web@tableX[,y]
    y.ind  <- web@dictionary[,y]
    y.vals <- `dimnames<-`(web@tableY[,y.ind], NULL)

    #straight to web!!
    disjointed[[t]]@estimates[web@dictionary[,y],y]  <- R_disj(y.vals, qr.mm, i.disj)

    emergent  [[t]]@estimates[web@dictionary[,y],y]  <- R_emerg(y.vals, qr.mm, i.disj)
  }

}

#Add F and P statistics
modelfit   <- lapply(modelfit,   get_PF, d = web@dictionary)
disjointed <- lapply(disjointed, get_PF, d = web@dictionary)
emergent   <- lapply(emergent,   get_PF, d = web@dictionary)


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
#'
dfmat <- function(x.assign, x.int, all.assign, x.fct, n){
  df0   <- colSums(!sapply(x.assign, function(x) index.self.high(x, all.assign, x.fct)))
  df1   <- c(sum(index.self.high(x.assign[1], all.assign, x.fct)),
             sapply(x.int, function(x) sum(all.assign %in% x)))
  df2 = n - (df0 + df1)
  dfr = df2/df1
  return(rbind(df1,df2, dfr))
}

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

R_disj_emerg <- function(y, x, mm, web, x.vars, i.disj){
  #adjust the input model.matrix by multiplying the relevant columns by x
  qr.mm  <- mm
  qr.mm[,x.vars] <- mm[,x.vars] * web@tableX[,y]
  y.ind  <- web@dictionary[,y]
  y.vals <- `dimnames<-`(web@tableY[,y.ind], NULL)

  Rsq_d <- R_disj(y.vals, qr.mm, i.disj)

  Rsq_e <- R_emerg(y.vals, qr.mm, i.disj)

  return(cbind(Rsq_d,
               Rsq_e))

}

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


  R_emerg <- function(y.vals, qr.mm, i.disj){
    #Emergent
    e.ord   <- order(rowSums(qr.mm[,i.disj[,3], drop = FALSE]))
    e.mm    <- cbind(1, diff(qr.mm[e.ord,-1]))

    e.y     <- diff(y.vals[e.ord,])

    #use diff to look at var
    e.mm.0  <- e.mm[,i.disj[,1]]
    e.mm.1  <- e.mm[,i.disj[,2]]

    #Cycle through dropping interactions with x, including higher order interactions.
    RSS_i0 <- SS(fast.qr.resid(y = e.y, x = e.mm.0))

    #Cycle through returning interactions with x, but not higher order interactions.
    RSS_i1 <- SS(fast.qr.resid(y = e.y, x = e.mm.1))

    return(1 - (RSS_i1 / RSS_i0))

  }



R_full <- function(y, mm, web, x.fct, x.vars){
  #adjust the input model.matrix by multiplying the relevant columns by x
  qr.mm  <- mm
  qr.mm[,x.vars] <- mm[,x.vars] * web@tableX[,y]
  y.ind  <- web@dictionary[,y]
  y.vals <- `dimnames<-`(web@tableY[,y.ind], NULL)

  #Return H1 SSR
  SS(fast.qr.resid(y = y.vals, x = qr.mm))

}

index.self.high <- function(x, all.assign, x.fct){
  all.assign %in% which(apply(
    x.fct * x.fct[,x] == x.fct[,x],
    MARGIN =  2, all))
}

index.high <- function(x, all.assign, x.fct){

  all.assign %in% which(
    `[[<-`(apply(
    x.fct * x.fct[,x] == x.fct[,x],
    MARGIN =  2, all), subscript = x, FALSE)
    )
}

index.disj <- function(x, all.assign, x.fct){
  #first is full null,
  #second is parameter to investigate returned
  i0 <- !all.assign %in% which(apply(
    x.fct * x.fct[,x] == x.fct[,x],
    MARGIN =  2, all))
  ix <- all.assign %in% x
  cbind(i0, i1 = i0 | ix, ix)
}

#' Calculate an association network
#' @description This is the main workspider function in the anansi package. It manages the individual functionalities of anansi, including correlation analysis, correlation by group and differential correlation.
#' @param web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.
#' @param method Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor()} are "spearman" and "kendall".
#' @param metadata A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
#' @param groups A vector of the column names of categorical values in the metadata object to control which groups should be assessed for simple correlations. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.
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
anansi2 = function(web, method = "pearson", groups = NULL,  metadata = NULL, formula = ~1,
                  adjust.method = "BH", modeltype = "lm", resampling = F, locality = F,
                  reff = NULL, verbose = T, diff_cor = T, ignore_dictionary = F){

  #Ensure appropriate options for type III ANOVAs
  contr.opt <- options(contrasts = c("contr.sum","contr.poly"))

  if(is.vector(metadata)){
    metadata = data.frame(metadata = metadata)
  }


  #If there is a formula that is not empty:
  if(!identical(all.vars(formula), character(0))){
    unique_factors = intersect(
      c(reff, all.vars(formula)),
      names(which(sapply(metadata,
                         function(x) is.character(x) || is.factor(x) || is.logical(x))))
    )

    if(!is.null(groups)){
      unique_factors = groups
    }

    metadata_factors = metadata[,unique_factors]
    groups = do.call(paste, c(data.frame(metadata_factors), sep = "_"))

    if(length(c(groups)) == 0){
      groups <- NULL
    }

  }

  if(!is.null(reff)){
    reff = metadata[,reff]
  }


  if(ignore_dictionary){
    if(verbose){print("Dictionary will be ignored. Running all vs all associations.")}
    #set dictionary to all TRUE
    web@dictionary <- web@dictionary == web@dictionary
    if(is(web, "argonansiWeb")){web@strat_dict <- web@strat_dict == web@strat_dict}
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
  outYarn = new("anansiYarn", input = new("anansiInput", web = web, groups = groups, formula = formula, reff = reff))

  if(verbose){print("Running annotation-based correlations")}

  #initialize cor_output list object
  output = new("anansiOutput")

  output@cor_results     = call_groupwise(web     = web, method  = method,
                                          groups  = groups, verbose = verbose)

  if(diff_cor){
    if(verbose){print("Fitting models for differential correlation testing")
      print(paste("Model type:", modeltype, sep = ""))}
    output@model_results = unlist(anansiDiffCor2(web = web, formula = formula,
                                                metadata = metadata, reff = reff,
                                                modeltype = modeltype, verbose = verbose))
  }
  outYarn@output = output

  #FDR
  outYarn <- anansiAdjustP(x = outYarn, method = adjust.method, resampling = resampling, locality = locality, verbose = verbose)

  #Restore default options
  on.exit(options(contr.opt))
  return(outYarn)
}

