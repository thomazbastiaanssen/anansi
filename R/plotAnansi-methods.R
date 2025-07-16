#' Dissociation plot
#'
#' `plotAnansi` generates an association plot from the output of
#' [anansi()] in the table format. It provides a convenient way to
#' visually assess relevant results from the anansi analysis, either in the
#' form of a dotplot or a graph.
#'
#' @param x a `data.frame` object output of [anansi()] in
#'   the table format.
#'
#' @param layout `Character scalar`. Specifies the plot layout to generate. It
#'   must be one of `c("dotplot, graph)`. (Default: `dotplot`)
#'
#' @param association.type `Character scalar`. Specifies the type of
#'   association to show in the plot. One of `"disjointed"`,
#'   `"emergent"` and `"full"`. (Default: `NULL`)
#'
#' @param model.var `Character scalar`. Specifies the name of a variable
#'   in the anansi model. It is relevant only when `association.type` is
#'   `"disjointed"` or `"emergent"`. (Default: `NULL`)
#'
#' @param group `Character scalar`. Selects one of the groups included in the
#'   anansi model. It is relevant only when `layout` is `graph`.
#'   (Default: `All`)
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
#'    should be filled (Default: `"group"`)
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
#' @param show.cor `Logical scalar`. Whether correlation edges should be
#'   labelled with correlation coefficients when `layout` is `graph`.
#'   (Default: `FALSE`)
#' @param ... additional parameters
#'
#' @details
#' `plotAnansi` provides a standardised method to visualise the results
#' of anansi by means of a differential association plot. The input for this
#' function should be generated from [anansi()] or
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
#' library(ggraph)
#'
#' web <- randomWeb(n_samples = 100)
#' mae <- as(web, "MultiAssayExperiment")
#'
#' # Perform anansi analysis
#' out <- weaveWeb(mae,
#'     tableY = "y", tableX = "x"
#' ) |> anansi(formula = ~group_ab)
#'
#' # Select significant interactions
#' out <- out[out$full_p.values < 0.05, ]
#'
#' # Visualise disjointed associations filled by group
#' plotAnansi(out,
#'     association.type = "disjointed",
#'     model.var = "group_ab",
#'     signif.threshold = 0.05,
#'     fill_by = "group"
#' )
#'
#' # Visualise full associations filled by group
#' plotAnansi(out,
#'     association.type = "full",
#'     signif.threshold = 0.05,
#'     fill_by = "group"
#' )
#'
#' # Visualise full associations as graph
#' plotAnansi(out,
#'     layout = "graph",
#'     association.type = "full",
#'     signif.threshold = 0.05,
#'     show.cor = TRUE
#' )
#'
#' # Visualise disjointed associations as graph
#' plotAnansi(out,
#'     layout = "graph",
#'     association.type = "disjointed",
#'     model.var = "group_ab",
#'     signif.threshold = 0.05
#' )
#'
#' @seealso [anansi()]
#' @name plotAnansi
#'
NULL

#' @name plotAnansi
#' @rdname plotAnansi
#' @export
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes theme guides labs geom_vline geom_point
#'     scale_x_continuous scale_y_discrete scale_alpha_manual theme_bw
#' @importFrom ggforce facet_col
#' @importFrom stats setNames
#' @importFrom S4Vectors isEmpty
method(plotAnansi, class_data.frame) <- function(
        x,
        layout = "dotplot",
        association.type = NULL,
        model.var = NULL,
        group = "All",
        signif.threshold = NULL,
        colour_by = NULL,
        color_by = colour_by,
        fill_by = "group",
        size_by = NULL,
        shape_by = NULL,
        y_position = "right",
        x_lab = "cor",
        y_lab = "",
        show.cor = FALSE
    ) {
        # Create list of Booleans whether args are defined
        defined_args <- lapply(
            list(
                association = association.type,
                model.var = model.var,
                signif = signif.threshold
            ),
            function(x) !is.null(x)
        )
        # Check association.type
        if (defined_args[["association"]]) {
            match.arg(
                association.type,
                choices = c("disjointed", "emergent", "full")
            )
        }
        # Check model.var
        if (defined_args[["model.var"]]) {
            match.arg(model.var, names(attr(x, "model_terms")))
        }
        # Check that group is valid for graph layout
        if (layout == "graph") {
            match.arg(group, attr(x, "group_terms")$All)
        }
        # Check association.type and model.var
        if (
            defined_args[["association"]] &&
                !defined_args[["model.var"]] &&
                association.type %in% c("disjointed", "emergent")
        ) {
            stop(
                "'model.var' must specify a variable of the anansi model ",
                "when 'association type' is set to ",
                association.type,
                call. = FALSE
            )
        }
        if (
            defined_args[["association"]] &&
                defined_args[["model.var"]] &&
                association.type == "full"
        ) {
            model.var <- NULL
            warning(
                "'model.var' is ignored when 'association type' ",
                "is set to ",
                association.type,
                call. = FALSE
            )
            model.var <- NULL
        }
        # Derive p-value column from association.type and model.var
        pval <-
            paste0(c(association.type, model.var, "p.values"), collapse = "_")
        # Check x
        if (isEmpty(x)) {
            stop("'x' is an empty data.frame", call. = FALSE)
        }
        if (!all(c("feature_X", "feature_Y") %in% colnames(x))) {
            stop(
                "'x' must be the output of 'anansi' in the table format ",
                "and must contain columns 'feature_X' ,'feature_Y'",
                call. = FALSE
            )
        }
        if (!any(grepl(pval, names(x)))) {
            stop("Could not find p-values in 'x'.", call. = FALSE)
        }
        # Convert anansi wide to long format
        x <- .wide2long(x)
        # Select group for graph layout
        if (layout == "graph") {
            x <- x[x$group == group, ]
        }
        # Update colour_by if color_by is defined
        if (!is.null(color_by) && is.null(colour_by)) {
            colour_by <- color_by
        }
        # Check aesthetics
        defined_args <- c(
            defined_args,
            mapply(
                .check_aes,
                aes_name = c("colour_by", "fill_by", "size_by", "shape_by"),
                aes_var = list(colour_by, fill_by, size_by, shape_by),
                MoreArgs = list(x = x),
                SIMPLIFY = FALSE
            )
        )
        # Check signif.threshold
        if (
            defined_args[["signif"]] &&
                (!is.numeric(signif.threshold) ||
                    signif.threshold < 0 ||
                    signif.threshold > 1)
        ) {
            stop(
                "'signif.threshold' must be a number between 0 and 1",
                call. = FALSE
            )
        }
        if (!defined_args[["association"]] && defined_args[["signif"]]) {
            warning(
                "'signif.threshold' is ignored when ",
                "'association type' is not defined",
                call. = FALSE
            )
        }
        # Check y_position
        match.arg(y_position, choices = c("left", "right"))
        # Check show.cor
        if (!is.logical(show.cor)) {
            stop("'show.cor' must be either TRUE or FALSE.", call. = FALSE)
        }
        # Assemble plot data
        pData <- data.frame(
            x = x[["r.values"]],
            y = x[["feature_X"]],
            facet = x[["feature_Y"]],
            colour = if (defined_args[["colour_by"]]) x[[colour_by]] else NA,
            fill = if (defined_args[["fill_by"]]) x[[fill_by]] else NA,
            size = if (defined_args[["size_by"]]) x[[size_by]] else NA,
            shape = if (defined_args[["shape_by"]]) x[[shape_by]] else NA,
            alpha = if (defined_args[["signif"]]) {
                factor(x[[pval]] < signif.threshold, levels = c(TRUE, FALSE))
            } else {
                NA
            }
        )
        if (layout == "dotplot") {
            # Generate dotplot
            p <- .create_dotplot(
                pData,
                defined_args,
                association.type,
                signif.threshold,
                colour_by,
                fill_by,
                shape_by,
                size_by,
                y_position,
                x_lab,
                y_lab
            )
        } else if (layout == "graph") {
            # Generate graph plot
            p <- .create_graphplot(
                pData,
                defined_args,
                association.type,
                signif.threshold,
                show.cor
            )
        }
        return(p)
}

################################ HELP FUNCTIONS ################################
# Convert anansi wide to long format
#' @description
#' Base R pivot longer for group terms of anansi output object
#' @param x data.frame, anansi() output.
#' @return a pivoted table
#' @noRd
.wide2long <- function(x) {
    mt <- attr(x, "model_terms")
    gt <- attr(x, "group_terms")
    groups <- gt$All
    # To dodge partial matches, require front.
    gr_regex <- paste0("^", groups, "_")

    gterms <- gsub("All_", "", colnames(x)[grepl("^All_", colnames(x))])
    l <- lapply(
        X = gr_regex,
        FUN = function(y) {
            `colnames<-`(
                x[, grepl(x = colnames(x), y), drop = FALSE],
                gterms
            )
        }
    )
    d <- do.call(rbind.data.frame, l)

    f <- `row.names<-.data.frame`(
        x[,
            -unlist(
                lapply(
                    X = gr_regex,
                    FUN = function(y) grep(x = colnames(x), y)
                ),
                FALSE,
                FALSE
            ),
            drop = FALSE
        ],
        NULL
    )
    # If "All" is in group column, make it the last level.
    groups <- factor(
        rep(groups, each = NROW(x)),
        levels = c(
            groups[groups != "All"],
            groups[groups == "All"]
        )
    )
    x <- cbind(f, group = groups, d)
    # Restore terms
    x <- `attr<-`(x, "model_terms", mt)
    x <- `attr<-`(x, "group_terms", gt)
    x
}

# Check aesthetics
.check_aes <- function(x, aes_name, aes_var) {
    # Check if aesthetic is defined
    aes_defined <- !is.null(aes_var)
    # Rise exception if aesthetic is not character or not in x
    if (aes_defined && !(aes_var %in% colnames(x) && is.character(aes_var))) {
        stop(
            "'",
            aes_name,
            "' must be a character string specifying the",
            " name of a 'groups' term used in the original anansi call",
            call. = FALSE
        )
    }
    return(aes_defined)
}
# Create graph plot
#' @importFrom patchwork wrap_plots plot_layout
.create_graphplot <- function(
    pData,
    defined_args,
    association.type,
    signif.threshold,
    show.cor
) {
  requireNamespace("ggraph")
    graph_list <- lapply(
        unique(pData$facet),
        .plot_facet_graph,
        pData = pData,
        defined_args = defined_args,
        association.type = association.type,
        signif.threshold = signif.threshold,
        show.cor = show.cor
    )

    p <- wrap_plots(graph_list) +
      plot_layout(guides = "collect")

    return(p)
}

#' @importFrom tidygraph tbl_graph
#' @importFrom ggraph ggraph geom_edge_link scale_edge_colour_gradient2
#'   geom_node_point geom_node_text theme_graph
.plot_facet_graph <- function(
    facet,
    pData,
    defined_args,
    association.type,
    signif.threshold,
    show.cor
) {
  requireNamespace("ggraph")
    # Select group
    pData <- pData[pData$facet == facet, ]
    # Retrieve node data
    node_data <- data.frame(
        node_key = c(facet, pData$y),
        alpha = factor(c(TRUE, as.character(pData$alpha)),
            levels = c(TRUE, FALSE))
    )
    # Retrieve edge data
    edge_data <- data.frame(
        from = facet,
        to = pData$y,
        cor = pData$x,
        label = if (show.cor) round(pData$x, 2) else NA
    )
    # Combine node and edge data into graph
    graph <- tbl_graph(nodes = node_data, edges = edge_data)
    # Visualise graph
    p <- ggraph(graph, layout = 'kk') +
        geom_edge_link(aes(colour = .data$cor, label = .data$label),
            label_size = 3) +
        scale_edge_colour_gradient2(low = "blue", high = "red", limits = c(-1, 1))
    # Add nodes
    if (defined_args[["signif"]]) {
        p <- p + geom_node_point(aes(alpha = .data$alpha))
    } else {
        p <- p + geom_node_point()
    }
    # Add node text
    p <- p + geom_node_text(aes(label = .data$node_key), repel = TRUE)
    # Add significance legend
    if (defined_args[["association"]] && defined_args[["signif"]]) {
        p <- p +
            scale_alpha_manual(
                values = c("TRUE" = 1, "FALSE" = 1 / 3),
                paste(association.type, "association\np <", signif.threshold)
            )
    }
    # Apply graph theme
    p <- p + theme_graph()
    # Remove legend when all associations are (non-)significant
    if (length(unique(node_data$alpha)) != 2L) {
        p <- p +
            guides(alpha = "none")
    }
    return(p)
}
# Create dotplot
.create_dotplot <- function(
    pData,
    defined_args,
    association.type,
    signif.threshold,
    colour_by,
    fill_by,
    shape_by,
    size_by,
    y_position,
    x_lab,
    y_lab
) {
    # Create base plot
    p <- ggplot(data = pData) +
        aes(
            x = .data$x,
            y = .data$y,
            colour = .data$colour,
            fill = .data$fill,
            shape = .data$shape,
            size = .data$size,
            alpha = .data$alpha
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "red")
    # Set point size, shape and border colour if not defined
    point_args <- list()
    if (!defined_args[["size_by"]]) {
        point_args["size"] <- 3
    }
    if (!defined_args[["shape_by"]]) {
        point_args["shape"] <- 21
    }
    if (!defined_args[["colour_by"]]) {
        point_args["colour"] <- "black"
    }
    # Add points and facets
    p <- p +
        do.call(geom_point, point_args) +
        facet_col(~ .data$facet, space = "free", scales = "free_y") +
        scale_x_continuous(
            limits = c(-1, 1),
            n.breaks = 11,
            expand = c(0, 0)
        ) +
        scale_y_discrete(limits = rev, position = y_position)
    # Add significance legend
    if (defined_args[["association"]] && defined_args[["signif"]]) {
        p <- p +
            scale_alpha_manual(
                values = c("TRUE" = 1, "FALSE" = 1 / 3),
                paste(association.type, "association\np <", signif.threshold)
            )
    }
    # Add labels
    p <- p +
        theme_bw() +
        labs(
            x = x_lab,
            y = y_lab,
            fill = fill_by,
            colour = colour_by,
            shape = shape_by,
            size = size_by
        )
    # Remove legend if aesthetics and significance are not defined
    if (
        !any(unlist(
            defined_args[c("colour_by", "fill_by", "size_by", "shape_by")]
        )) &&
            !defined_args[["signif"]]
    ) {
        p <- p + theme(legend.position = "none")
    }
    # Remove legend for undefined aesthetics
    guide_names <- c("colour", "fill", "alpha")[
        !unlist(defined_args[c("colour_by", "fill_by", "signif")])
    ]
    guide_args <- setNames(
        as.list(rep("none", length(guide_names))),
        guide_names
    )
    p <- p + do.call(guides, guide_args)

    return(p)
}
