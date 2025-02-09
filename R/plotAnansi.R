#' Dissociation plot
#' 
#' \code{plotAnansi} generates a standard dissociation plot from the output of
#' \code{\link{getAnansi}} in the long format. It provides a convenient way to
#' visually assess relevant results from the anansi analysis.
#' 
#' @param df a \code{data.frame} object output of \code{\link{getAnansi}} in
#'   the long format.
#'   
#' @param x \code{Character scalar}. Specifies the type of coefficient to
#' show in the plot. One of \code{"r.values"}, \code{"r.squared"} and
#'   \code{"q.values"}. (Default: \code{"r.values"})
#' 
#' @param association.type \code{Character scalar}. Specifies the type of
#' association to show in the plot. One of \code{"disjointed"},
#'   \code{"emergent"} and \code{"full"}. (Default: \code{"disjointed"})
#' 
#' @param model.var \code{Character scalar}. Specifies the name of a variable
#' in the anansi model. It is relevant only when \code{association.type} is
#'   \code{"disjointed"} or \code{"emergent"}. (Default: \code{NULL})
#'   
#' @param signif.threshold \code{Numeric scalar}. Specifies the significance
#'   threshold to show in the plot. (Default: \code{NULL})
#'
#' @param colour_by \code{Character scalar}. Specifies the name of a column in
#'   \code{df} by which points should be coloured. (Default: \code{NULL})
#'
#' @param color_by \code{Character scalar}. Alias to \code{colour_by}.
#' 
#' @param fill_by \code{Character scalar}. Specifies the name of a column in
#'   \code{df} by which points should be filled (Default: \code{NULL})
#' 
#' @param size_by \code{Character scalar}. Specifies the name of a column in
#'   \code{df} by which points should be sized (Default: \code{NULL})
#' 
#' @param shape_by \code{Character scalar}. Specifies the name of a column in
#'   \code{df} by which points should be shaped (Default: \code{NULL})
#' 
#' @param x_lab \code{Character scalar}. Specifies the label of the x axis.
#'   (Default: \code{"cor"})
#' 
#' @param y_lab \code{Character scalar}. Specifies the label of the y axis.
#'   (Default: \code{""})
#' 
#' @param y_position \code{Character scalar}. Specifies the position of the y
#'   labels. It should be either \code{"left"} or \code{"right"}.
#'   (Default: \code{"right"})
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
#' # Import libraries
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
    function(df, x = "r.values", association.type = "disjointed", model.var = NULL,
    colour_by = NULL, color_by = colour_by, fill_by = NULL, size_by = NULL,
    shape_by = NULL, signif.threshold = NULL, y_position = "right",
    x_lab = "cor", y_lab = ""){
      
    match.arg(x, choices = c("r.values", "r.squared", "q.values"))
    match.arg(association.type, choices = c("disjointed", "emergent", "full"))
    
    if( association.type %in% c("disjointed", "emergent") && is.null(model.var) ){
        stop("'model.var' must specify a variable of the anansi model when ",
            "'association type' is set to ", association.type, call. = FALSE)
    }
    if( association.type == "full" && !is.null(model.var) ){
        warning("'model.var' is ingored when 'association type' is set to ",
            association.type, call. = FALSE)
    }
      
    pval <- paste0(c("model", association.type, model.var, "p.values"),
        collapse = "_")
    x <- paste0(c("model", association.type, model.var, x), collapse = "_")
    
    # Check df
    if( !all(c(x, pval, "feature_X", "feature_Y") %in% colnames(df)) ){
        stop("'df' must be the output of 'getAnansi' in the long format and ",
        "must contain columns '", x, "' ,'", pval, "', 'feature_X', 'feature_Y'",
        call. = FALSE)
    }
    # Check y_position
    match.arg(y_position, choices = c("left", "right"))
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
    pData <- data.frame(x = df[[x]], y = df[["feature_X"]],
        colour = if (colour_defined) df[[colour_by]] else NA,
        fill = if (fill_defined) df[[fill_by]] else NA,
        size = if (size_defined) df[[size_by]] else NA,
        shape = if (shape_defined) df[[shape_by]] else NA,
        alpha = if (signif_defined)
            factor(df[[pval]] < signif.threshold,
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
            paste(association.type, "association\np <", signif.threshold))
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