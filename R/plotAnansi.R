#' Dissociation plot
#'
#' \code{plotAnansi} generates a standard dissociation plot from the output of
#' \code{\link{getAnansi}} in the table format. It provides a convenient way to
#' visually assess relevant results from the anansi analysis.
#'
#' @param x a \code{data.frame} object output of \code{\link{getAnansi}} in
#'   the table format.
#'  
#' @param association.type \code{Character scalar}. Specifies the type of
#' association to show in the plot. One of \code{"disjointed"},
#'   \code{"emergent"} and \code{"full"}. (Default: \code{NULL})
#'
#' @param model.var \code{Character scalar}. Specifies the name of a variable
#' in the anansi model. It is relevant only when \code{association.type} is
#'   \code{"disjointed"} or \code{"emergent"}. (Default: \code{NULL})
#'
#' @param signif.threshold \code{Numeric scalar}. Specifies the threshold to
#'   mark the significance of \code{association.type}. (Default: \code{NULL})
#'
#' @param colour_by \code{Character scalar}. Specifies the name of a column in
#'   \code{x} by which points should be coloured. (Default: \code{NULL})
#'
#' @param color_by \code{Character scalar}. Alias to \code{colour_by}.
#'
#' @param fill_by \code{Character scalar}. Specifies the name of a column in
#'   \code{x} by which points should be filled (Default: \code{NULL})
#'
#' @param size_by \code{Character scalar}. Specifies the name of a column in
#'   \code{x} by which points should be sized (Default: \code{NULL})
#'
#' @param shape_by \code{Character scalar}. Specifies the name of a column in
#'   \code{x} by which points should be shaped (Default: \code{NULL})
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
#' of anansi by means of a differential association plot. The input for this
#' function should be generated from \code{\link{getAnansi}} or
#' \code{\link{anansi}}, with \code{return.format = "table"}
#'
#' @return
#' A ggplot2 object.
#'
#' @examples
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
#' out <- out[out$full_q.values < 0.1, ]
#'
#' # Visualise disjointed associations coloured by group
#' plotAnansi(out,
#'            association.type = "disjointed",
#'            model.var = "Legend",
#'            signif.threshold = 0.05,
#'            fill_by = "group")
#'            
#' # Visualise full associations coloured by group
#' plotAnansi(out,
#'            association.type = "full",
#'            signif.threshold = 0.05,
#'            fill_by = "group")
#'            
#' @seealso
#' \code{\link{getAnansi}}
#' \code{\link{anansi}}
#'
#' @name plotAnansi
#'
NULL

#' @rdname plotAnansi
#' @export
setGeneric("plotAnansi", signature = c("x"),
    function(x, ...) standardGeneric("plotAnansi")
)

#' @rdname plotAnansi
#' @export
#' @importFrom ggplot2 ggplot theme guides labs geom_vline geom_point
#'   scale_x_continuous scale_y_discrete scale_alpha_manual theme_bw
#' @importFrom ggforce facet_col
#' @importFrom stats setNames
#' @importFrom S4Vectors isEmpty
setMethod("plotAnansi", signature = c(x = "data.frame"),
    function(x, association.type = NULL, model.var = NULL,
    signif.threshold = NULL, colour_by = NULL, color_by = colour_by,
    fill_by = NULL, size_by = NULL, shape_by = NULL, y_position = "right",
    x_lab = "cor", y_lab = ""){
    # Create list of Booleans whether args are defined
    defined_args <- lapply(list(association = association.type,
        model.var = model.var, signif = signif.threshold),
        function(x) !is.null(x))
    # Check association.type
    if( defined_args[["association"]] ){
        match.arg(association.type, choices = c("disjointed", "emergent", "full"))
    }
    # Check model.var
    if( defined_args[["model.var"]] ){
        match.arg(model.var, attr(x, "model_terms"))
    }
    # Check association.type and model.var
    if( defined_args[["association"]] && !defined_args[["model.var"]] &&
        association.type %in% c("disjointed", "emergent") ){
        stop("'model.var' must specify a variable of the anansi model when ",
            "'association type' is set to ", association.type, call. = FALSE)
    }
    if( defined_args[["association"]] && defined_args[["model.var"]] &&
        association.type == "full" ){
        model.var <- NULL
        warning("'model.var' is ignored when 'association type' is set to ",
            association.type, call. = FALSE)
        model.var <- NULL
    }
    # Derive p-value column from association.type and model.var
    pval <- paste0(c(association.type, model.var, "p.values"), collapse = "_")
    # Check x
    if( isEmpty(x) ){
        stop("'x' is an empty data.frame", call. = FALSE)
    }
    if( !all(c("feature_X", "feature_Y") %in% colnames(x)) ||
        !any(grepl(pval, names(x)))){
        stop("'x' must be the output of 'anansi' in the table format and must ",
        "contain columns 'feature_X' ,'feature_Y', 'r.values' and '", pval, "'",
        call. = FALSE)
    }
    # Convert anansi wide to long format
    x <- .wide2long(x)
    # Update colour_by if color_by is defined
    if (!is.null(color_by) && is.null(colour_by)) {
        colour_by <- color_by
    }
    # Check aesthetics
    defined_args <- c(defined_args,
        mapply(.check_aes,
            aes_name = c("colour_by", "fill_by", "size_by", "shape_by"),
            aes_var = list(colour_by, fill_by, size_by, shape_by),
            MoreArgs = list(x = x), SIMPLIFY = FALSE))
    # Check signif.threshold
    if( defined_args[["signif"]] && (!is.numeric(signif.threshold)
        || signif.threshold < 0 || signif.threshold > 1) ){
        stop("'signif.threshold' must be a number between 0 and 1", call. = FALSE)
    }
    if( !defined_args[["association"]] && defined_args[["signif"]] ){
        warning("'signif.threshold' is ignored when 'association type' is not",
            " defined", call. = FALSE)
    }
    # Check y_position
    match.arg(y_position, choices = c("left", "right"))
    # Assemble plot data
    pData <- data.frame(x = x[["r.values"]], y = x[["feature_X"]],
        colour = if (defined_args[["colour_by"]]) x[[colour_by]] else NA,
        fill = if (defined_args[["fill_by"]]) x[[fill_by]] else NA,
        size = if (defined_args[["size_by"]]) x[[size_by]] else NA,
        shape = if (defined_args[["shape_by"]]) x[[shape_by]] else NA,
        alpha = if (defined_args[["signif"]])
            factor(x[[pval]] < signif.threshold,
            levels = c(TRUE, FALSE)) else NA,
        facet = x[["feature_Y"]]
    )
    # Generate dotplot
    p <- .create_dotplot(pData, defined_args, association.type,
        signif.threshold, colour_by, fill_by, shape_by, size_by, y_position,
        x_lab, y_lab)
    return(p)
})
################################ HELP FUNCTIONS ################################
# Convert anansi wide to long format
#' @importFrom tidyr pivot_longer pivot_wider separate_wider_regex
.wide2long <- function(x){
    # Fetch group terms
    group_terms <- attr(x, "group_terms")
    # Create regex pattern to match group terms
    patterns <- paste0("^(?:", paste(group_terms, collapse = "|"), ")")
    # Convert anansi wide to long format
    x_long <- x |>
        pivot_longer(matches(patterns)) |>
        separate_wider_regex(name,
            c(group = patterns, "_", cor_param = ".\\.values$")) |>
        pivot_wider(names_from = cor_param, values_from = value)
    return(x_long)
}
# Check aesthetics
.check_aes <- function(x, aes_name, aes_var){
    # Check if aesthetic is defined
    aes_defined <- !is.null(aes_var)
    # Rise exception if aesthetic is not character or not in x
    if( aes_defined && !(aes_var %in% colnames(x) && is.character(aes_var)) ){
        stop("'", aes_name, "' must be a character string specifying the",
            " name of a column in 'x'", call. = FALSE)
    }
    return(aes_defined)
}
# Create dotplot
.create_dotplot <- function(pData, defined_args, association.type,
        signif.threshold, colour_by, fill_by, shape_by, size_by, y_position,
        x_lab, y_lab){
    # Create base plot
    p <- ggplot(data = pData, aes(x = .data$x, y = .data$y,
            colour = .data$colour, fill = .data$fill, shape = .data$shape,
            size = .data$size, alpha = .data$alpha)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "red")
    # Set point size and shape if not defined
    point_args <- list()
    if( !defined_args[["size_by"]] ){
        point_args["size"] <- 3
    }
    if( !defined_args[["shape_by"]] ){
        point_args["shape"] <- 21
    }
    # Add points and facets
    p <- p + do.call(geom_point, point_args) +
        facet_col(~ .data$facet, space = "free", scales = "free_y") +
        scale_x_continuous(limits = c(-1, 1), n.breaks = 11, expand = c(0, 0)) +
        scale_y_discrete(limits = rev, position = y_position)
    # Add significance legend
    if( defined_args[["association"]] && defined_args[["signif"]] ){
        p <- p + scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1/3),
            paste(association.type, "association\np <", signif.threshold))
    }
    # Add labels
    p <- p + theme_bw() +
        labs(x = x_lab, y = y_lab, fill = fill_by, colour = colour_by,
            shape = shape_by, size = size_by)
    # Remove legend if aesthetics and significance are not defined
    if( !any(unlist(defined_args[c("colour_by", "fill_by", "size_by", "shape_by")])) &&
        !defined_args[["signif"]] ){
        p <- p + theme(legend.position = "none")
    }
    # Remove legend for undefined aesthetics
    guide_names <- c("colour", "fill", "alpha")[
        !unlist(defined_args[c("colour_by", "fill_by", "signif")])
    ]
    guide_args <- setNames(as.list(rep("none", length(guide_names))), guide_names)
    p <- p + do.call(guides, guide_args)
    return(p)
}