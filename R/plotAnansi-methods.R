#' Dissociation plot
#'
#' `plotAnansi` generates a standard dissociation plot from the output of
#' [getAnansi()] in the table format. It provides a convenient way to
#' visually assess relevant results from the anansi analysis.
#'
#' @param x a `data.frame` object output of [getAnansi()] in
#'   the table format.
#'
#' @param association.type `Character scalar`. Specifies the type of
#' association to show in the plot. One of `"disjointed"`,
#'   `"emergent"` and `"full"`. (Default: `NULL`)
#'
#' @param model.var `Character scalar`. Specifies the name of a variable
#' in the anansi model. It is relevant only when `association.type` is
#'   `"disjointed"` or `"emergent"`. (Default: `NULL`)
#'
#' @param signif.threshold `Numeric scalar`. Specifies the threshold to
#'   mark the significance of `association.type`. (Default: `NULL`)
#'
#' @param colour_by `Character scalar`. Specifies one of the `groups`
#'    terms used in the original `anansi` call, `x` by which points
#'    should be coloured. (Default: `NULL`)
#'
#' @param color_by `Character scalar`. Alias to `colour_by`.
#'
#' @param fill_by `Character scalar`. Specifies one of the `groups`
#'    terms used in the original `anansi` call, `x` by which points
#'    should be filled (Default: `NULL`)
#'
#' @param size_by `Character scalar`. Specifies one of the `groups`
#'    terms used in the original `anansi` call, `x` by which points
#'    should be sized. (Default: `NULL`)
#'
#' @param shape_by `Character scalar`. Specifies one of the `groups`
#'    terms used in the original `anansi` call, `x` by which points
#'    should be shaped. (Default: `NULL`)
#'
#' @param x_lab `Character scalar`. Specifies the label of the x axis.
#'   (Default: `"cor"`)
#'
#' @param y_lab `Character scalar`. Specifies the label of the y axis.
#'   (Default: `""`)
#'
#' @param y_position `Character scalar`. Specifies the position of the y
#'   labels. It should be either `"left"` or `"right"`.
#'   (Default: `"right"`)
#'
#' @param ... additional arguments
#'
#' @details
#' `plotAnansi` provides a standardised method to visualise the results
#' of anansi by means of a differential association plot. The input for this
#' function should be generated from [getAnansi()] or
#' [anansi()], with `return.format = "table"`
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
#' web <- randomWeb(n_samples = 100)
#' mae <- as(web, "MultiAssayExperiment")
#'
#' # Perform anansi analysis
#' out <- getAnansi(mae,
#'   tableY = "y", tableX = "x",
#'   formula = ~ cat_ab
#' )
#'
#' # Select significant interactions
#' out <- out[out$full_p.values < 0.05, ]
#'
#' # Visualise disjointed associations filled by group
#' plotAnansi(out,
#'            association.type = "disjointed",
#'            model.var = "cat_ab",
#'            signif.threshold = 0.05,
#'            fill_by = "group")
#'
#' # Visualise full associations filled by category ('cat_ab')
#' plotAnansi(out,
#'            association.type = "full",
#'            signif.threshold = 0.05,
#'            fill_by = "cat_ab")
#'
#' @seealso
#' [getAnansi()]
#' [anansi()]
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
#' @importFrom ggplot2 ggplot aes .data theme guides labs geom_vline geom_point
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
        match.arg(model.var, names(attr(x, "model_terms")))
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
#' @importFrom dplyr across mutate
#' @importFrom tidyr pivot_longer pivot_wider separate_wider_regex replace_na
#' @importFrom tidyselect all_of matches
.wide2long <- function(x){
    # Create regex pattern to match group terms
    patterns <- lapply(attr(x, "group_terms"),
        function(term) paste0("(?:", paste(term, collapse = "|"), ")_?"))
    group_cols <- c("group", names(patterns[-1]))
    # Convert anansi wide to long format
    x_long <- x |>
        pivot_longer(matches(paste0("^", patterns[["All"]]))) |>
        separate_wider_regex(name,
            c(group = patterns[["All"]], cor_param = ".\\.values$")) |>
        pivot_wider(names_from = cor_param, values_from = value) |>
    # Separate the group column into the original group terms
        separate_wider_regex(group, unlist(patterns[-1]),
            too_few = "align_start", cols_remove = FALSE) |>
        mutate(across(all_of(group_cols),
            ~ gsub("_$", "", replace_na(., "All"))))
    return(x_long)
}
# Check aesthetics
.check_aes <- function(x, aes_name, aes_var){
    # Check if aesthetic is defined
    aes_defined <- !is.null(aes_var)
    # Rise exception if aesthetic is not character or not in x
    if( aes_defined && !(aes_var %in% colnames(x) && is.character(aes_var)) ){
        stop("'", aes_name, "' must be a character string specifying the",
            " name of a 'groups' term used in the original anansi call",
            call. = FALSE)
    }
    return(aes_defined)
}
# Create dotplot
.create_dotplot <- function(pData, defined_args, association.type,
        signif.threshold, colour_by, fill_by, shape_by, size_by, y_position,
        x_lab, y_lab){
    # Create base plot
    p <- ggplot(data = pData) +
      aes(x = .data$x, y = .data$y, colour = .data$colour,
          fill = .data$fill, shape = .data$shape,
          size = .data$size, alpha = .data$alpha) +
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
