#' Take anansi output and wrangle it all to a wide data.frame format.
#' @param anansi_output The output of the main anansi function.
#' @param prune Boolean, default is TRUE. Toggles whether to take out the non-canonical associations.
#' @param translate Boolean, default is FALSE Toggles whether to translate the names of the features in tableX and tableY to human readable names.
#' @param Y_translation data.frame, a lookup table with featureY names as the first column and human readable names as the second. See \code{cpd_translation} for an example file.
#' @param X_translation data.frame, a lookup table with featureX names as the first column and human readable names as the second. See \code{KO_translation} for an example file.
#' @return a wide format data.frame
#' @export
#'
spinToWide <- function(anansi_output, prune = TRUE, translate = FALSE, Y_translation = NULL, X_translation = NULL) {
  # in case of argonansi
  if (is(anansi_output@input@web, "argonansiWeb")) {
    anansi_output@output@model_results[["modelfit.full"]] <- argonaut_rep_to_strat(
      x = anansi_output@output@model_results[["modelfit.full"]],
      web = anansi_output@input@web
    )
  }
  # First flatten all types  of results (cor, model, etc) and create a list of individual wide data.frames.
  flat_list <- lapply(unlist(list(
    anansi_output@output@cor_results,
    anansi_output@output@model_results
  )), getAnansiResults)



  # merge all data.frames in the list while keeping the row order.
  wide_df <- Reduce(function(x, y) merge(x, y, sort = FALSE), flat_list)


  if (translate) {
    argonansi <- FALSE
    if (is(anansi_output@input@web, "argonansiWeb")) {
      argonansi <- TRUE
    }
    wide_df <- anansiTranslate(wide_df, Y_translation = Y_translation, X_translation = X_translation, argonansi = argonansi)
  }


  if (prune) {
    if (is(anansi_output@input@web, "argonansiWeb")) {
      # If true, remove all non-canonical interactions.
      wide_df <- wide_df[c(anansi_output@input@web@strat_dict), ]
    }
    if (!is(anansi_output@input@web, "argonansiWeb")) {
      # If true, remove all non-canonical interactions.
      wide_df <- wide_df[c(anansi_output@input@web@dictionary), ]
    }
  }

  # Sort the output by the first column
  wide_df <- wide_df[order(wide_df[, 1]), ]

  return(wide_df)
}


#' Take anansi output and wrangle it all to a long data.frame format.
#' @param anansi_output an \code{anansiYarn} object. The output of the main anansi function.
#' @param prune Boolean, default is TRUE. Toggles whether to take out the non-canonical associations.
#' @param translate Boolean, default is FALSE Toggles whether to translate the names of the features in tableX and tableY to human readable names.
#' @param Y_translation data.frame, a lookup table with featureY names as the first column and human readable names as the second. See \code{cpd_translation} for an example file.
#' @param X_translation data.frame, a lookup table with featureX names as the first column and human readable names as the second. See \code{KO_translation} for an example file.
#' @return a long format data.frame intended to be compatible with \code{ggplot2}
#' @export
#'
spinToLong <- function(anansi_output, prune = TRUE, translate = FALSE, Y_translation = NULL, X_translation = NULL) {
  # Figure out how many groups were present in the correlations.
  n_cors <- length(anansi_output@output@cor_results)

  # Flatten all types  of the correlations and create a list of individual wide data.frames.
  flat_cor_list <- lapply(anansi_output@output@cor_results, getAnansiResults, format = "long")

  # Merge all individual correlation results into a single long data.frame
  long_out <- Reduce(rbind, flat_cor_list)
  colnames(long_out)[3] <- "r.values"

  if (length(anansi_output@output@model_results) > 0) {
    # Make a flat list of the model results in wide format
    # Catch argonansi exception
    if (is(anansi_output@input@web, "argonansiWeb")) {
      anansi_output@output@model_results[["modelfit.full"]] <- argonaut_rep_to_strat(
        x = anansi_output@output@model_results[["modelfit.full"]],
        web = anansi_output@input@web
      )
    }

    flat_model_list <- lapply(anansi_output@output@model_results, getAnansiResults, format = "wide")

    # Merge all model results in a single wide data.frame
    wide_model_df <- Reduce(function(x, y) merge(x, y, sort = FALSE), flat_model_list)

    # Stack the model results so that they have the same number of rows as the correlation results
    # wide_model_df   = do.call(rbind, replicate(n_cors, wide_model_df, simplify = FALSE))

    # Combine the results into a single data.frame
    long_out <- cbind(long_out, wide_model_df[, -c(1, 2)])
  }

  if (translate) {
    argonansi <- FALSE
    if (is(anansi_output@input@web, "argonansiWeb")) {
      argonansi <- TRUE
    }
    long_out <- anansiTranslate(long_out, Y_translation = Y_translation, X_translation = X_translation, argonansi = argonansi)
  }

  if (prune) {
    if (is(anansi_output@input@web, "argonansiWeb")) {
      # If true, remove all non-canonical interactions.
      long_out <- long_out[c(anansi_output@input@web@strat_dict), ]
    }
    if (!is(anansi_output@input@web, "argonansiWeb")) {
      # If true, remove all non-canonical interactions.
      long_out <- long_out[c(anansi_output@input@web@dictionary), ]
    }
  }

  return(long_out)
}

#' Extract information from an anansiTale object and parse it into a neat table
#' @param tale An \code{anansiTale} object
#' @param format either \code{wide} or \code{long}. Default is \code{wide}. Toggles the output format.
#' @return A wide format data.frame. The first two columns are the features from tableY and tableX, respectively
#' @seealso \code{\link{spinToWide}}\cr \code{\link{spinToLong}}
#'
getAnansiResults <- function(tale, format = "wide") {
  if (!format %in% c("wide", "long")) {
    stop("format should be either `wide` or `long`.")
  }
  out_df <-
    cbind(
      expand.grid(
        feature_Y = row.names(tale@estimates),
        feature_X = colnames(tale@estimates), stringsAsFactors = FALSE
      ),
      estimates = c(tale@estimates),
      p.values = c(tale@p.values),
      q.values = c(tale@q.values)
    )

  if (format == "wide") {
    names(out_df)[3] <- tale@type
    names(out_df)[-c(1,2)] <- paste(tale@subject, names(out_df)[-c(1,2)], sep = "_")
  }
  if (format == "long") {
    out_df$type <- paste(tale@subject, sep = "_")
  }

  return(out_df)
}


#' Translate tableY and tableX columns column to human readable names
#' @description helper function for the spinToLong and \code{spinToWide} functions.
#' @param x a table with a feature_Y and feature_X column.
#' @param Y_translation data.frame, a lookup table with featureY names as the first column and human readable names as the second. See \code{cpd_translation} for an example file.
#' @param X_translation data.frame, a lookup table with featureX names as the first column and human readable names as the second. See \code{KO_translation} for an example file.
#' @param argonansi A boolean. Toggles whether we are dealing with stratified data.
#' @return an expanded table with a column for the human readable names for the features in both tableY and tableX.
#'
anansiTranslate <- function(x, Y_translation = Y_translation, X_translation = X_translation, argonansi = FALSE) {
  l_X <- length(unique(x$feature_X))
  l_Y <- rep(length(unique(x$feature_Y)), l_X)
  if (argonansi) {
    x$feature_X_stratified <- x$feature_X
    x$feature_X <- gsub(x = x$feature_X_stratified, pattern = "\\..*", replacement = "")
    x$subtype_X <- gsub(x = x$feature_X_stratified, pattern = ".*\\.", replacement = "")
    l_Y <- rle(paste(x$feature_X))$lengths[seq_len(length(unique(x$feature_X)))]
  }

  relevant_Y <- Y_translation[Y_translation[, 1] %in% x$feature_Y, ]
  x$feature_Y <- rep(
    paste(relevant_Y[, 1],
      gsub(";.*", "", relevant_Y[, 2]),
      sep = ": "
    ),
    times = l_X
  )

  relevant_X <- X_translation[X_translation[, 1] %in% x$feature_X, ]
  x$feature_X <- rep(
    paste(relevant_X[, 1],
      relevant_X[, 2],
      sep = ": "
    ),
    l_Y
  )
  return(x)
}

#' Take anansi output and wrangle it to a list of plottable objects.
#' @param anansiYarn The output of the main anansi function.
#' @param target A boolean matrix. Determines which associations should be prepared for plotting. Default is the dictionary.
#' @param translate Boolean, default is FALSE Toggles whether to translate the names of the features in tableX and tableY to human readable names.
#' @param Y_translation data.frame, a lookup table with featureY names as the first column and human readable names as the second. See \code{cpd_translation} for an example file.
#' @param X_translation data.frame, a lookup table with featureX names as the first column and human readable names as the second. See \code{KO_translation} for an example file.
#' @return a list of ready to plot data.frames and their names.
#' @export
#' @examples
#' # Starting off from the example in ?anansi
#'
#' data(dictionary)
#' data(FMT_data)
#'
#' # Clean and prepare the example data.
#' # In the example dataset, the metabolites are already cleaned.
#'
#' KOs <- floor(FMT_KOs)
#' KOs <- apply(KOs, c(1, 2), function(x) as.numeric(as.character(x)))
#' KOs <- KOs[apply(KOs == 0, 1, sum) <= (ncol(KOs) * 0.90), ]
#' KOs.exp <- clr_c(KOs[row.names(KOs) %in% sort(unique(unlist(anansi_dic))), ])
#'
#' t1 <- t(FMT_metab)
#' t2 <- t(KOs.exp)
#'
#' # Run anansi pipeline.
#'
#' web <- weaveWebFromTables(
#'   tableY = t1,
#'   tableX = t2,
#'   dictionary = anansi_dic
#' )
#'
#' anansi_out <- anansi(
#'   web = web,
#'   formula = ~Legend,
#'   groups = FMT_metadata$Legend,
#'   metadata = FMT_metadata,
#'   adjust.method = "BH",
#'   verbose = TRUE
#' )
#'
#' # Let's look at all canonical interactions that also have a sufficiently well-fitting model:
#'
#' outPlots <- spinToPlots(anansi_out)
#'
#' # load ggplot2 and patchwork for plotting
#'
#' library(ggplot2)
#' library(patchwork)
#'
#' plotted <- lapply(outPlots, FUN = function(p) {
#'   # Main ggplot call
#'   ggplot(data = p$data, aes(x = X, y = Y, fill = groups)) +
#'
#'     # Define geoms:
#'     geom_point(shape = 21) +
#'     geom_smooth(method = "lm") +
#'     theme_bw() +
#'
#'     # Improve annotation:
#'     scale_fill_manual(values = c(
#'       "Young yFMT" = "#2166ac",
#'       "Aged oFMT" = "#b2182b",
#'       "Aged yFMT" = "#ef8a62"
#'     )) +
#'     ylab(p$name[1]) +
#'     xlab(p$name[2]) +
#'     ggtitle(paste(p$name[1], "vs", p$name[2]))
#' })
#'
#' # Call patchwork to unify and arrange the first 6 plots
#'
#' wrap_plots(plotted[1:6]) + plot_layout(guides = "collect")
#'
spinToPlots <- function(anansiYarn, target = NULL, translate = FALSE, Y_translation = NULL, X_translation = NULL) {

  if(is.null(target)) {target = yarn.dic.logical(anansiYarn)} else {
    target = yarn.dic.logical(anansiYarn)& target}

   pairs_of_interest <- which(target, arr.ind = TRUE, useNames = FALSE)

  out_list <- apply(X = pairs_of_interest, MARGIN = 1, FUN = gather_plot, anansiYarn = anansiYarn, simplify = FALSE)

  names(out_list) <- paste(
    colnames(anansiYarn@input@web@tableY)[pairs_of_interest[, 1]],
    "vs",
    colnames(anansiYarn@input@web@tableX)[pairs_of_interest[, 2]]
  )

  # Reorganize order of out_list as to reflect order of the other spinToX functions.
  out_list <- out_list[order(names(out_list))]

  if (translate) {
    out_list <- lapply(out_list, FUN = function(y) {
      y$name <- unlist(anansiTranslate(
        x = tiny_df(y),
        Y_translation = Y_translation,
        X_translation = X_translation
      ))
      return(y)
    })
  }

  return(out_list)
}

#' Take anansi output and wrangle it to a list of plottable objects. Meant to be apply'd by \code{spinToPlots()}.
#' @param pair_of_interest A vector of length 2. Denotes the column position of the tableY and tableX data to be plotted.
#' @param anansiYarn The output of the main anansi function.
#' @return a list of ready to plot data.frames and their names
#'
gather_plot <- function(pair_of_interest, anansiYarn) {
  df_out <- data.frame(
    Y = anansiYarn@input@web@tableY[, pair_of_interest[1]],
    X = anansiYarn@input@web@tableX[, pair_of_interest[2]],
    groups = anansiYarn@input@groups
  )
  return(list(
    data = df_out,
    name = c(
      colnames(anansiYarn@input@web@tableY)[pair_of_interest[1]],
      colnames(anansiYarn@input@web@tableX)[pair_of_interest[2]]
    )
  ))
}

#' Helper function to translate feature names in spinToPlot()
#' @param x a vector of two names, to be wrangled to a 1x2 data.frame.
#' @noRd
#'
tiny_df <- function(x) {
  return(data.frame(matrix(x$name, nrow = 1, dimnames = list("", c("feature_Y", "feature_X")))))
}
