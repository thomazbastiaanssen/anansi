#' Handle FDR methods for anansi
#'
anansiPAdjust <- function(p, dictionary, method, best.of.two){
  if(best.of.two){

    p = q

    return(q)
  }

  else{return(p.adjust(p[dictionary], method = method))}
}


#
#
# a1 = t(apply(X = p[d], MARGIN = c(1), FUN = function(x){p.adjust(p = x, method = "BH")}))
# a2 = (apply(X = p, MARGIN = c(2), FUN = function(x){p.adjust(p = x, method = "BH")}))
# a3 = (apply(X = p, MARGIN = c(1,2), FUN = function(x){p.adjust(p = x, method = "BH")}))
#
# p[,1]
#
# exp(mean(log(c(1, 10))))
#
# p.adjust(c(0.1, 0.5), method = "BH", n = 3)
#
#
#
# View(rbind(a1, a2))
# View(do.call(rbind, list(p, a1, a2, a3)))
#
# p = anansi_out@output@cor_results$`Aged yFMT`@p.values
# d = anansi_out@input@web@dictionary
#
# ?(p.adjust)
#
#
# data(dictionary)
# data(FMT_data)
#
# #Clean and prepare the example data.
# #In the example dataset, the metabolites are already cleaned.
#
# KOs   <- floor(FMT_KOs)
# KOs   <- apply(KOs,c(1,2),function(x) as.numeric(as.character(x)))
# KOs   <- KOs[apply(KOs == 0, 1, sum) <= (ncol(KOs) * 0.90), ]
#
# KOs   <- KOs[row.names(KOs) %in% sort(unique(unlist(anansi_dic))),]
#
# #CLR-transform.
#
# KOs.exp = clr_c(KOs)
#
# #Make sure that columns are features and rows are samples.
#
# t1 = t(FMT_metab)
# t2 = t(KOs.exp)
#
# #Run anansi pipeline.
#
# web        = weaveWebFromTables(tableY     = t1,
#                                 tableX     = t2,
#                                 dictionary = anansi_dic)
#
# anansi_out = anansi(web     = web,
#                     method  = "pearson",
#                     groups  = FMT_metadata$Legend,
#                     adjust.method = "BH",
#                     verbose = TRUE)
