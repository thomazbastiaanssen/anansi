#' Methods for AnansiWeb S7 container class
#' @name AnansiWeb-methods
#' @param x input, AnansiWeb object
#' @aliases show.anansi::AnansiWeb names.anansi::AnansiWeb
#' @aliases dimnames.anansi::AnansiWeb dim.anansi::AnansiWeb
#' @seealso \itemize{
#' \item [weaveWeb()]: for general use.
#' \item [AnansiWeb-pairwise]: for methods for pairwise operations
#' }
#' @returns The desired information from an AnansiWeb object
#' @examples
#' # Setup
#' web <- randomWeb(n_samp = 36)
#'
#' # Accessors
#' dimnames(web)
#' dim(web)
#' names(web)
#'
#' # Getters and setters:
#'
#' tableX(web)[1:5, 1:5]
#' tableY(web)[1:5, 1:5]
#' dictionary(web)
#' head(metadata(web))
#'
#' # Assign some random metadata
#' metadata(web) <- data.frame(
#'     id = row.names(tableY(web)),
#'     a = rnorm(36),
#'     b = sample(c("a", "b"), 36, TRUE),
#'     row.names = "id"
#' )
#'
#' # Coerce to list
#' weblist <- as.list(web)
#'
#' # Coerce to Data.frame
#' webdf <- as.data.frame(web)
#'
#' # Coerce to MultiAssayExperiment
#' mae <- asMAE(web)
#'
#' # Coerce to TreeSummarizedExperiment
#' tse <- asTSE(web)
NULL

#' @importFrom methods show
#' @importMethodsFrom methods show
#' @export
#'
S7::method(show, AnansiWeb) <- function(object) {
    cat(
        paste(class(object), collapse = " "),
        " with ",
        NROW(object@tableX),
        " observations:\n    tableY: ",
        names(object)[1],
        " (",
        NROW(object),
        " features)\n    tableX: ",
        names(object)[2],
        " (",
        NCOL(object),
        " features)\n",
        sep = ""
    )
    cat("Use @ to access: tableX, tableY, dictionary, metadata.")
    invisible(NULL)
}

#' @export
#'
S7::method(print, AnansiWeb) <- function(x, ...) show(x)


#' @export
#'
S7::method(dimnames, AnansiWeb) <- function(x) dimnames(x@dictionary)

#' @export
#'
S7::method(dim, AnansiWeb) <- function(x) dim(x@dictionary)

S7::method(names, AnansiWeb) <- function(x) names(dimnames(x@dictionary))

#' @export
#'
S7::method(tableX, AnansiWeb) <- function(x) S7::prop(x, "tableX")

#' @export
#'
S7::method(tableY, AnansiWeb) <- function(x) S7::prop(x, "tableY")

#' @export
#'
S7::method(dictionary, AnansiWeb) <- function(x) S7::prop(x, "dictionary")

#' Get metadata.
#' @name metadata
#' @details
#' Compatible with S4Vectors generic.
#' @rdname metadata
#' @importMethodsFrom S4Vectors metadata
#' @aliases metadata,anansi::AnansiWeb-method
#' @param x input object
#' @param ... additional arguments
#' @returns `metadata` slot of x
#' @examples
#' x <- randomWeb(10)
#' metadata(x)
#' @export
#'
S7::method(metadata, AnansiWeb) <- function(x, ...) S7::prop(x, "metadata")

#' Set metadata.
#' @name metadata<-
#' @aliases metadata.set
#' @details
#' Compatible with S4Vectors generic.
#' @importMethodsFrom S4Vectors metadata<-
#' @aliases metadata<-,anansi::AnansiWeb-method
#' @param x input object
#' @param ... additional arguments
#' @param value replacement value, coerced to data.frame
#' @returns `x` with modified `metadata` slot.
#' @examples
#' x <- randomWeb(10)
#' metadata(x) <- cbind(
#'     metadata(x),
#'     new_groups = c("A", "B")
#' )
#' @export
#'
S7::method(`metadata<-`, AnansiWeb) <- function(x, ..., value) {
    x@metadata <- as.data.frame(value)
    x
}

################################################################################
################################################################################

#' For S4Vectors::metadata compatibility
#' @noRd
metadata <- function(x) S7::prop(x, "metadata")

#' For S4Vectors::metadata<- compatibility
#' @noRd
`metadata<-` <- function(x, value) {
    xx@metadata <- as.data.frame(value)
    x
}

#' Is this a data.frame with exactly two columns that are named?
#' @noRd
validWeb <- function(x) {
    y_names <- identical(rownames(x), colnames(x@tableY))
    x_names <- identical(colnames(x), colnames(x@tableX))
    s_names <- identical(rownames(x@tableY), rownames(x@tableX))
    meta_dim <- any(
        NROW(x@metadata) == NROW(x@tableY),
        prod(dim(x@metadata)) <= 1
    )
    if (!y_names) {
        message("colnames(tableY), rownames(dictionary) not identical.")
    }
    if (!x_names) {
        message("colnames(tableX), colnames(dictionary) not identical.")
    }
    if (!s_names) {
        message("rownames(tableX), rownames(tableX) not identical.")
    }
    if (!meta_dim) {
        message("NROW(metadata) does not equal rows of tableY, tableX.")
    }

    return(all(y_names, y_names, s_names, meta_dim))
}
