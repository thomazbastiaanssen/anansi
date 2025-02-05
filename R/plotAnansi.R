#' Dissociation plot
#' 
#' @param df
#'
#' @param colour_by
#'
#' @param color_by
#' 
#' @param fill_by
#' 
#' @param size_by
#' 
#' @param shape_by
#' 
#' @param signif.threshold
#' 
#' @param x_lab
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
#'   scale_x_continuous scale_y_discrete scale_alpha_manual element_blank
#' @importFrom ggforce facet_col
setMethod("plotAnansi", signature = c(df = "data.frame"),
    function(df, colour_by = NULL, color_by = colour_by, fill_by = NULL,
    size_by = NULL, shape_by = NULL, signif.threshold = NULL, x_lab = "cor") {
    
    if( !all(c("r.values", "feature_X", "feature_Y") %in% colnames(df)) ){
        stop("'df' must be the output of 'getAnansi' in the long format and ",
        "must contain columns 'r.values', 'feature_X', 'feature_Y' and ",
        "'model_disjointed_Legend_p.values'", call. = FALSE)
    }
      
    if (!is.null(color_by) && is.null(colour_by)) {
        colour_by <- color_by
    }
    
    colour_defined <- .check_aes(df, colour_by, "colour_by")
    fill_defined <- .check_aes(df, fill_by, "fill_by")
    size_defined <- .check_aes(df, size_by, "size_by")
    shape_defined <- .check_aes(df, shape_by, "shape_by")
    
    signif_defined <- !is.null(signif.threshold)
    if( signif_defined && !is.numeric(signif.threshold) ){
        stop("'signif.threshold' must be a number between 0 and 1", call. = FALSE)
    }
    
    pData <- data.frame(x = df[["r.values"]], y = df[["feature_X"]],
        colour = NA, fill = NA, size = NA, shape = NA, alpha = NA,
        facet = df[["feature_Y"]])
  
    if( colour_defined ){
        pData["colour"] <- df[[colour_by]]
    }
    if( fill_defined ){
        pData["fill"] <- df[[fill_by]]
    }
    if( size_defined ){
        pData["size"] <- df[[size_by]]
    }
    if( shape_defined ){
        pData["shape"] <- df[[shape_by]]
    }
    
    if( signif_defined ){
        pData["alpha"] <- df["model_disjointed_Legend_p.values"] < signif.threshold
        pData["alpha"] <- factor(pData[["alpha"]], levels = c(TRUE, FALSE))
    }
    
    p <- ggplot(data = pData, aes(x = .data$x, y = .data$y,
            colour = .data$colour, fill = .data$fill, shape = .data$shape,
            size = .data$size, alpha = .data$alpha)) + 
        geom_vline(xintercept = 0, linetype = "dashed", colour = "red")
    
    if( !size_defined && !shape_defined ){
        p <- p + geom_point(shape = 21, size = 3)
    } else if( !size_defined && shape_defined ){
        p <- p + geom_point(size = 3)
    } else if( size_defined && !shape_defined ){
        p <- p + geom_point(shape = 21)
    } else if( size_defined && shape_defined ){
        p <- p + geom_point()
    }
    
    p <- p +
        facet_col(~ .data$facet, space = "free", scales = "free_y") +
        scale_x_continuous(limits = c(-1, 1), n.breaks = 11, expand = c(0, 0)) +
        scale_y_discrete(limits = rev, position = "right")
    
    if( signif_defined ){
        p <- p + scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1/3),
            paste("Disjointed association\np <", signif.threshold))
    }
    
    p <- p + theme_bw() +
        labs(x = x_lab, fill = fill_by, colour = colour_by,
            shape = shape_by, size = size_by) +
        theme(axis.title.y = element_blank())
    
    if( !(colour_defined || fill_defined || size_defined || shape_defined) && !signif_defined ){
        p <- p + theme(legend.position = "none")
    }
    if( !colour_defined ){
        p <- p + guides(colour = "none")
    }
    if( !fill_defined ){
        p <- p + guides(fill = "none")
    }
    if( !signif_defined ){
        p <- p + guides(alpha = "none")
    }
    
    return(p)
})


.check_aes <- function(df, aes_var, aes_name) {
  
    aes_defined <- !is.null(aes_var)
    
    if( aes_defined && !(aes_var %in% colnames(df) && is.character(aes_var)) ){
        stop("'", aes_name, "' must be a character string specifiying the",
            " name of a column in 'df'", call. = FALSE)
    }

    return(aes_defined)
}
  
