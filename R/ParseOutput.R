#' Extract information from an anansiTale object and parse it into a neat table
#' @param tale An \code{anansiTale} object
#' @param dic A dictionary.
#' @return A wide format data.frame with summary statistics by feature pair.
#'
frame.tale <- function(tale, dic){
  switch(tale@type,
         "r.values"  = frame.tale.cor(tale, dic),
         "r.squared" = frame.tale.ols(tale, dic))
}

#' @noRd
#'
frame.tale.cor <- function(tale, dic){
  out.df <- data.frame(
    r.values = tale@estimates[dic],
    t.values = tale@t.values[dic],
    p.values = tale@p.values[dic]
    #q.values =  tale@q.values[dic]
    )
  colnames(out.df) <- paste(tale@subject, colnames(out.df), sep = "_")
  return(out.df)
}

#' @noRd
#'
frame.tale.ols <- function(tale, dic){
  out.df <- data.frame(
    r.squared = tale@estimates[dic],
    f.values =  tale@f.values[dic],
    p.values =  tale@p.values[dic]
    #q.values =  tale@q.values[dic]
    )
  colnames(out.df) <- paste(tale@subject, colnames(out.df), sep = "_")
  return(out.df)
}

#' Shape list of `anansiTale` results into a data.frame.
#' @noRd
#'
result.df <- function(out.list, dic) {
  feature_labs <- expand.grid(
         feature_Y = row.names(dic),
         feature_X = colnames(dic), stringsAsFactors = FALSE
     )[dic,]

  df.list <- c(feature_labs, lapply(out.list, frame.tale, dic))
  do.call(what = "cbind.data.frame", args = df.list, quote = TRUE)
}

# ## Translate tableY and tableX columns column to human readable names
# ## @description helper function for the spinToLong and \code{spinToWide} functions.
# ## @param x a table with a feature_Y and feature_X column.
# ## @param Y_translation data.frame, a lookup table with featureY names as the first column and human readable names as the second. See \code{cpd_translation} for an example file.
# ## @param X_translation data.frame, a lookup table with featureX names as the first column and human readable names as the second. See \code{KO_translation} for an example file.
# ## @param argonansi A boolean. Toggles whether we are dealing with stratified data.
# ## @return an expanded table with a column for the human readable names for the features in both tableY and tableX.
# ##
# anansiTranslate <- function(x, Y_translation = Y_translation, X_translation = X_translation, argonansi = FALSE) {
#   l_X <- length(unique(x$feature_X))
#   l_Y <- rep(length(unique(x$feature_Y)), l_X)
#   if (argonansi) {
#     x$feature_X_stratified <- x$feature_X
#     x$feature_X <- gsub(x = x$feature_X_stratified, pattern = "\\..*", replacement = "")
#     x$subtype_X <- gsub(x = x$feature_X_stratified, pattern = ".*\\.", replacement = "")
#     l_Y <- rle(paste(x$feature_X))$lengths[seq_len(length(unique(x$feature_X)))]
#   }
#
#   relevant_Y <- Y_translation[Y_translation[, 1] %in% x$feature_Y, ]
#   x$feature_Y <- rep(
#     paste(relevant_Y[, 1],
#           gsub(";.*", "", relevant_Y[, 2]),
#           sep = ": "
#     ),
#     times = l_X
#   )
#
#   relevant_X <- X_translation[X_translation[, 1] %in% x$feature_X, ]
#   x$feature_X <- rep(
#     paste(relevant_X[, 1],
#           relevant_X[, 2],
#           sep = ": "
#     ),
#     l_Y
#   )
#   return(x)
# }
#
# ## Take anansi output and wrangle it to a list of plottable objects.
# ## @param anansiYarn The output of the main anansi function.
# ## @param target A boolean matrix. Determines which associations should be prepared for plotting. Default is the dictionary.
# ## @param translate Boolean, default is FALSE Toggles whether to translate the names of the features in tableX and tableY to human readable names.
# ## @param Y_translation data.frame, a lookup table with featureY names as the first column and human readable names as the second. See \code{cpd_translation} for an example file.
# ## @param X_translation data.frame, a lookup table with featureX names as the first column and human readable names as the second. See \code{KO_translation} for an example file.
# ## @return a list of ready to plot data.frames and their names.
# ##
# ## #examples
# ## # Starting off from the example in ?anansi
# ##
# ## data(dictionary)
# ## data(FMT_data)
# ##
# ## # Clean and prepare the example data.
# ## # In the example dataset, the metabolites are already cleaned.
# ##
# ## KOs <- floor(FMT_KOs)
# ## KOs <- apply(KOs, c(1, 2), function(x) as.numeric(as.character(x)))
# ## KOs <- KOs[apply(KOs == 0, 1, sum) <= (ncol(KOs) * 0.90), ]
# ## KOs.exp <- clr_c(KOs[row.names(KOs) %in% sort(unique(unlist(anansi_dic))), ])
# ##
# ## t1 <- t(FMT_metab)
# ## t2 <- t(KOs.exp)
# ##
# ## # Run anansi pipeline.
# ##
# ## web <- weaveWebFromTables(
# ##   tableY = t1,
# ##   tableX = t2,
# ##   dictionary = anansi_dic
# ## )
# ##
# ## anansi_out <- anansi(
# ##   web = web,
# ##   formula = ~Legend,
# ##   groups = "Legend",
# ##   metadata = FMT_metadata,
# ##   adjust.method = "BH",
# ##   verbose = TRUE
# ## )
# ##
# ## # Let's look at all canonical interactions that also have a sufficiently well-fitting model:
# ##
# ## outPlots <- spinToPlots(anansi_out)
# ##
# ## # load ggplot2 and patchwork for plotting
# ##
# ## library(ggplot2)
# ## library(patchwork)
# ##
# ## plotted <- lapply(outPlots, FUN = function(p) {
# ##   # Main ggplot call
# ##   ggplot(data = p$data, aes(x = X, y = Y, fill = groups)) +
# ##
# ##     # Define geoms:
# ##     geom_point(shape = 21) +
# ##     geom_smooth(method = "lm") +
# ##     theme_bw() +
# ##
# ##     # Improve annotation:
# ##     scale_fill_manual(values = c(
# ##       "Young yFMT" = "#2166ac",
# ##       "Aged oFMT" = "#b2182b",
# ##       "Aged yFMT" = "#ef8a62"
# ##     )) +
# ##     ylab(p$name[1]) +
# ##     xlab(p$name[2]) +
# ##     ggtitle(paste(p$name[1], "vs", p$name[2]))
# ## })
# ##
# ## # Call patchwork to unify and arrange the first 6 plots
# ##
# ## wrap_plots(plotted[1:6]) + plot_layout(guides = "collect")
# ## @noRd
# spinToPlots <- function(anansiYarn, target = NULL, translate = FALSE, Y_translation = NULL, X_translation = NULL) {
#   if (is.null(target)) {
#     target <- yarn.dic.logical(anansiYarn)
#   } else {
#     target <- yarn.dic.logical(anansiYarn) & target
#   }
#
#   pairs_of_interest <- which(target, arr.ind = TRUE, useNames = FALSE)
#
#   out_list <- apply(X = pairs_of_interest, MARGIN = 1, FUN = gather_plot, anansiYarn = anansiYarn, simplify = FALSE)
#
#   names(out_list) <- paste(
#     colnames(anansiYarn@input@web@tableY)[pairs_of_interest[, 1]],
#     "vs",
#     colnames(anansiYarn@input@web@tableX)[pairs_of_interest[, 2]]
#   )
#
#   # Reorganize order of out_list as to reflect order of the other spinToX functions.
#   out_list <- out_list[order(names(out_list))]
#
#   if (translate) {
#     out_list <- lapply(out_list, FUN = function(y) {
#       y$name <- unlist(anansiTranslate(
#         x = tiny_df(y),
#         Y_translation = Y_translation,
#         X_translation = X_translation
#       ))
#       return(y)
#     })
#   }
#
#   return(out_list)
# }
#
# ## Take anansi output and wrangle it to a list of plottable objects. Meant to be apply'd by \code{spinToPlots()}.
# ## @param pair_of_interest A vector of length 2. Denotes the column position of the tableY and tableX data to be plotted.
# ## @param anansiYarn The output of the main anansi function.
# ## @return a list of ready to plot data.frames and their names
# ##
# gather_plot <- function(pair_of_interest, anansiYarn) {
#   df_out <- data.frame(
#     Y = anansiYarn@input@web@tableY[, pair_of_interest[1]],
#     X = anansiYarn@input@web@tableX[, pair_of_interest[2]],
#     groups = anansiYarn@input@groups
#   )
#   return(list(
#     data = df_out,
#     name = c(
#       colnames(anansiYarn@input@web@tableY)[pair_of_interest[1]],
#       colnames(anansiYarn@input@web@tableX)[pair_of_interest[2]]
#     )
#   ))
# }
#
# ## Helper function to translate feature names in spinToPlot()
# ## @param x a vector of two names, to be wrangled to a 1x2 data.frame.
# ## @noRd
# ##
# tiny_df <- function(x) {
#   return(data.frame(matrix(x$name, nrow = 1, dimnames = list("", c("feature_Y", "feature_X")))))
# }
