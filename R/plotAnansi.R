#' Dissociation plot
#' 
#' \code{plotAnansi} generates a standard dissociation plot from the output of
#' \code{\link{getAnansi}} in the long format. It provides a convenient way to
#' visually assess relevant results from the anansi analysis.
#' 
#' @param df a \code{data.frame} object output of \code{\link{getAnansi}} in
#'   the long format.
#'
#' @param colour_by \code{Character scalar}. Specifies the name of a column in
#' \code{df} by which points should be coloured. (Default: \code{NULL})
#'
#' @param color_by \code{Character scalar}. Alias to \code{colour_by}.
#' 
#' @param fill_by \code{Character scalar}. Specifies the name of a column in
#' \code{df} by which points should be filled (Default: \code{NULL})
#' 
#' @param size_by \code{Character scalar}. Specifies the name of a column in
#' \code{df} by which points should be sized (Default: \code{NULL})
#' 
#' @param shape_by \code{Character scalar}. Specifies the name of a column in
#' \code{df} by which points should be shaped (Default: \code{NULL})
#' 
#' @param signif.threshold \code{Numeric scalar}. Specifies the significance
#'   threshold to show in the plot. (Default: \code{NULL})
#' 
#' @param y_position \code{Character scalar}. Specifies the position of the y
#'   labels. It should be either \code{"left"} or \code{"right"}.
#'   (Default: \code{"right"})
#' 
#' @param x_lab \code{Character scalar}. Specifies the label of the x axis.
#'   (Default: \code{"cor"})
#' 
#' @param y_lab \code{Character scalar}. Specifies the label of the y axis.
#'   (Default: \code{""})
#' 
#' @details
#' \code{plotAnansi} provides a standardised method to visualise the results
#' of anansi by means of a dissociation plot. The input for this function should
#' be generated from \code{\link{getAnansi}} with \code{return.format = "long"}
#' or from \code{\link{anansi}} followed by \code{\link{spinToLong}}.
#'
#' @return
#' A ggplot2 object.
#'
#' @examples
#' 
#' #' # Import libraries
#' library(mia)
#' library(TreeSummarizedExperiment)
#' library(MultiAssayExperiment)
#'
#' # Load data
#' data("FMT_data", package = "anansi")
#' data("dictionary", package = "anansi")
#'
#' # Convert to (Tree)SummarizedExperiment objects
#' metab_se <- SummarizedExperiment(assays = SimpleList(conc = as.matrix(FMT_metab)))
#' KO_tse <- TreeSummarizedExperiment(assays = SimpleList(counts = as.matrix(FMT_KOs)))
#'
#' # Select functions that are represented in the dictionary
#' keep <- row.names(KO_tse) %in% sort(unique(unlist(anansi_dic)))
#' KO_tse <- KO_tse[keep, ]
#'
#' # Remove features with less than 10% prevalence
#' KO_tse <- subsetByPrevalent(KO_tse,
#'   assay.type = "counts",
#'   prevalence = 0.1
#' )
#'
#' # Perform a centered log-ratio transformation on the functional counts assay
#' KO_tse <- transformAssay(KO_tse,
#'   assay.type = "counts",
#'   method = "clr",
#'   pseudocount = TRUE
#' )
#'
#' # Prepare colData
#' coldata <- FMT_metadata
#' rownames(coldata) <- coldata$Sample_ID
#' coldata <- coldata[match(colnames(KO_tse), rownames(coldata)), ]
#'
#' # Combine experiments into MultiAssayExperiment object
#' mae <- MultiAssayExperiment(
#'   experiments = ExperimentList(metabolites = metab_se, functions = KO_tse),
#'   colData = coldata
#' )
#'
#' # Perform anansi analysis
#' out <- getAnansi(mae,
#'   experiment1 = "metabolites", experiment2 = "functions",
#'   assay.type1 = "conc", assay.type2 = "clr",
#'   formula = ~Legend, translate = TRUE,
#'   X_translation = KO_translation, Y_translation = cpd_translation
#' )
#'
#' # Select significant interactions
#' out <- out[out$model_full_q.values < 0.1, ]
#'
#' # Generate dissociation plot
#' plotAnansi(out, fill_by = "type", signif.threshold = 0.05)
#' 
#' @seealso
#' \code{\link{getAnansi}}
#' \code{\link{anansi}}
#' \code{\link{spinToLong}}
#' 
#' @name plotAnansi
#' 
NULL

#' @rdname plotAnansi
#' @export
setGeneric("plotAnansi", signature = c("df"),
    function(df, ...) standardGeneric("plotAnansi")
)

#' @rdname plotAnansi
#' @export
#' @importFrom ggplot2 ggplot theme guides labs geom_vline geom_point
#'   scale_x_continuous scale_y_discrete scale_alpha_manual theme_bw
#' @importFrom ggforce facet_col
#' @importFrom stats setNames
setMethod("plotAnansi", signature = c(df = "data.frame"),
    function(df, colour_by = NULL, color_by = colour_by, fill_by = NULL,
    size_by = NULL, shape_by = NULL, signif.threshold = NULL,
    y_position = "right", x_lab = "cor", y_lab = "") {
    # Check df
    if( !all(c("r.values", "feature_X", "feature_Y", "model_disjointed_Legend_p.values") %in% colnames(df)) ){
        stop("'df' must be the output of 'getAnansi' in the long format and ",
        "must contain columns 'r.values', 'feature_X', 'feature_Y' and ",
        "'model_disjointed_Legend_p.values'", call. = FALSE)
    }
    # Update colour_by if color_by is defined
    if (!is.null(color_by) && is.null(colour_by)) {
        colour_by <- color_by
    }
    # Check aesthetics
    colour_defined <- .check_aes(df, colour_by, "colour_by")
    fill_defined <- .check_aes(df, fill_by, "fill_by")
    size_defined <- .check_aes(df, size_by, "size_by")
    shape_defined <- .check_aes(df, shape_by, "shape_by")
    # Check signif.threshold
    signif_defined <- !is.null(signif.threshold)
    if( signif_defined && !is.numeric(signif.threshold) ){
        stop("'signif.threshold' must be a number between 0 and 1", call. = FALSE)
    }
    # Assemple plot data
    pData <- data.frame(x = df[["r.values"]], y = df[["feature_X"]],
        colour = if (colour_defined) df[[colour_by]] else NA,
        fill = if (fill_defined) df[[fill_by]] else NA,
        size = if (size_defined) df[[size_by]] else NA,
        shape = if (shape_defined) df[[shape_by]] else NA,
        alpha = if (signif_defined)
            factor(df["model_disjointed_Legend_p.values"] < signif.threshold,
            levels = c(TRUE, FALSE)) else NA,
        facet = df[["feature_Y"]]
    )
    # Create base plot
    p <- ggplot(data = pData, aes(x = .data$x, y = .data$y,
            colour = .data$colour, fill = .data$fill, shape = .data$shape,
            size = .data$size, alpha = .data$alpha)) + 
        geom_vline(xintercept = 0, linetype = "dashed", colour = "red")
    # Set point size and shape if not defined
    point_args <- list()
    if( !size_defined ){
        point_args["size"] <- 3
    }
    if( !shape_defined ){
        point_args["shape"] <- 21
    }
    # Add points and facets
    p <- p + do.call(geom_point, point_args) +
        facet_col(~ .data$facet, space = "free", scales = "free_y") +
        scale_x_continuous(limits = c(-1, 1), n.breaks = 11, expand = c(0, 0)) +
        scale_y_discrete(limits = rev, position = y_position)
    # Add significance legend
    if( signif_defined ){
        p <- p + scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1/3),
            paste("Disjointed association\np <", signif.threshold))
    }
    # Add labels
    p <- p + theme_bw() +
        labs(x = x_lab, y = y_lab, fill = fill_by, colour = colour_by,
            shape = shape_by, size = size_by)
    # Remove legend if aesthetics and significance are not defined
    if( !(colour_defined || fill_defined || size_defined || shape_defined) && !signif_defined ){
        p <- p + theme(legend.position = "none")
    }
    # Remove legend for undefined aesthetics
    guide_names <- c("colour", "fill", "alpha")[!c(colour_defined, fill_defined, signif_defined)]
    guide_args <- setNames(as.list(rep("none", length(guide_names))), guide_names)
    p <- p + do.call(guides, guide_args)
    
    return(p)
})

################################ HELP FUNCTIONS ################################
# Check aesthetics
.check_aes <- function(df, aes_var, aes_name) {
    # Check if aesthetic is defined
    aes_defined <- !is.null(aes_var)
    # Rise exception if aesthetic is not character or not in df
    if( aes_defined && !(aes_var %in% colnames(df) && is.character(aes_var)) ){
        stop("'", aes_name, "' must be a character string specifiying the",
            " name of a column in 'df'", call. = FALSE)
    }
    return(aes_defined)
}