#' Handle FDR methods for anansi. Can also be used on anansi output to recalculate FDR.
#' @param x an AnansiYarn object, the main output from anansi().
#' @param method The p-value adjustment method. See ?p.adjust.
#' @param resampling A boolean. Toggles the resampling of p-values to help reduce the influence of dependence of p-values. Will take more time on large datasets.
#' @param locality A boolean. Toggles whether to prefer sampling from p-values from a comparison that shares an x or y feature. In a nutshell, considers the local p-value landscape when more important when correcting for FDR.
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @importFrom stats p.adjust
#' @importFrom methods slot slot<-
#' @export
#'
anansiAdjustP <- function(x, method = "BH", resampling = F, locality = T, verbose = T){
  #first check if the settings make sense
  stopifnot("Input needs to be an anansiYarn object, the output of the anansi() function" = class(x) == "anansiYarn")

  if(!resampling & locality){
    locality = F
    if(verbose){print("Locality is only available for resampling procedure. Setting locality to False")}
  }

  #verbalise the plan if verbose
  if(verbose)
  {if(method %in% c("BH", "fdr"))
    {print("Adjusting p-values using Benjamini & Hochberg's procedure.")}
    else
    {print("Adjusting p-values...")}

    if(resampling & locality){print("Using locally resampled distributions. ")}
    else if(resampling & !locality){print("Using total resampled distributions.")}
    else if(!resampling){print("Using theoretical distribution.")}

    }

  #extract dictionary for convenience
  dictionary = x@input@web@dictionary
  if(is(x@input@web, "argonansiWeb")){
    dictionary = x@input@web@strat_dict
  }

  #Take inventory of where to find p-values.
  pval_df = data.frame(model = c(
    rep("cor_results",   length(x@output@cor_results  )),
    rep("model_results", length(x@output@model_results))),
                  group = c(
                    names(x@output@cor_results),
                    names(x@output@model_results)))

  if(is(x@input@web, "argonansiWeb")){
    pval_df = pval_df[pval_df$group != "modelfit.full",]
    p = slot(slot(x@output, 'model_results')[['modelfit.full']], 'p.values')

    #Resampling approach
    if(resampling){
      #Extract coordinates to cycle over relevant p-values using apply
      coord <- which(x@input@web@dictionary, arr.ind = T)

      #adjust p-values directly into the relevant slot
      slot(slot(x@output, 'model_results')[['modelfit.full']], 'q.values')[x@input@web@dictionary] <-
        apply(coord,
              MARGIN = 1,
              FUN = function(y){
                compute_FDR(
                  coord      = y,
                  p          = p,
                  dictionary = x@input@web@dictionary,
                  locality   = locality,
                  method     = method)})
    }
    else{
      #adjust p-values directly into the relevant slot
      slot(slot(x@output, 'model_results')[['modelfit.full']], 'q.values')[x@input@web@dictionary] <-
        p.adjust(p[x@input@web@dictionary], method = method)
    }

  }

  #for each of those sources of p-values, do:
  for(i in 1:nrow(pval_df)){
    #Extract p-values
    p = slot(slot(x@output, pval_df[i,1])[[pval_df[i,2]]], 'p.values')

    #Resampling approach
    if(resampling){
      #Extract coordinates to cycle over relevant p-values using apply
      coord <- which(dictionary, arr.ind = T)

      #adjust p-values directly into the relevant slot
      slot(slot(x@output, pval_df[i,1])[[pval_df[i,2]]], 'q.values')[dictionary] <-
        apply(coord,
                     MARGIN = 1,
                     FUN = function(x){
                       compute_FDR(
                         coord      = x,
                         p          = p,
                         dictionary = dictionary,
                         locality   = locality,
                         method     = method)})
    }
    else{
      #adjust p-values directly into the relevant slot
      slot(slot(x@output, pval_df[i,1])[[pval_df[i,2]]], 'q.values')[dictionary] <-
      p.adjust(p[dictionary], method = method)
    }

  }
  return(x)

}


#' Determines which p-values to sample from
#' @param y the row-coordinate of the p-value of interest
#' @param x the column-coordinate of the p-value of interest
#' @param p A matrix containing the p-values from an anansiTale object.
#' @param dictionary And anansi Dictionary object.
#' @param locality A boolean. Toggles whether to prefer sampling from p-values from a comparison that shares an x or y feature. In a nutshell, considers the local p-value landscape when more important when correcting for FDR.
#'
determine_source <- function(y, x, p, dictionary, locality = T){
  if(!locality){return(unname(p[dictionary]))}
  else{
    return(unname(
      c(
      p[y,dictionary[y,]],
      p[dictionary[,x],x],
      median(p[dictionary])
      )
      )
    )
  }

}

#' Create a resampling matrix for FDR
#' @param coord a vector of two integers. the first contains the row and the second the columns coordinates of the p-value of interest.
#' @param p A matrix containing the p-values from an anansiTale object.
#' @param dictionary And anansi Dictionary object.
#' @param method The p-value adjustment method. See ?p.adjust.
#' @param locality A boolean. Toggles whether to prefer sampling from p-values from a comparison that shares an x or y feature. In a nutshell, considers the local p-value landscape when more important when correcting for FDR.
#'
compute_FDR <- function(coord, p, dictionary, locality, method){

  #Extract y (row) and x (column) coordinates of the relevant p-value.

  y = unlist(coord[1])
  x = unlist(coord[2])
  source = determine_source(x = x, y = y, p = p, dictionary = dictionary, locality = locality)

  #Generate resampling matrix, as to keep resampling consistent between p-values
  resam_mat = replicate(n = 1000,
                        sample(x       = source,
                               size    = length(p[dictionary]),
                               replace = T)
                        )

  resam_mat[,1] <- p[y, x]

  #Perform FDR on each row, only take the first position and return the median.
  q <- median(apply(X =  resam_mat,
                    MARGIN = 1,
                    FUN = function(x){p.adjust(x, method = method)})[1,]
              )

  return(q)
}
