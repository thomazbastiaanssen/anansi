#' AnansiWeb S7 container class
#' @name AnansiWeb
#' @description
#' `AnansiWeb` is an S7 class containing two feature tables as well as a
#' dictionary to link them. `AnansiWeb` is the main container that will
#' hold your input data throughout the `anansi` pipeline.
#'
#' Typical use of the `anansi` package will involve generating an `AnansiWeb`
#' object using the `weaveWeb()` function.
#'
#' The function `AnansiWeb()` constructs an `AnansiWeb` object from two
#' feature tables and an adjacency matrix.
#'
#' @param x input, AnansiWeb object
#' @seealso \itemize{
#' \item [weaveWeb()]: for general use.
#' \item [AnansiWeb-pairwise]: for methods for pairwise operations
#' }
#' @examples
#'
#' # Methods for AnansiWeb
#' dimnames(web)
#' dim(web)
#' names(web)
#'
#' web@tableX
#' web@tableY
#' web@dictionary
#'
#' # Assign some random metadata
#' web@metadata <- data.frame(
#'     id = row.names(web@tableY),
#'     a = rnorm(36),
#'     b = sample(c("a", "b"), 36, TRUE),
#'     row.names = "id"
#' )
#'
#' # coerce To list
#' weblist <- as.list(web)
#'
#' # Coerce to MultiAssayExperiment
#' asMAE(web)
NULL

#' @name show.AnansiWeb
#' @importFrom methods show
#' @importMethodsFrom methods show
#' @aliases show,anansi::AnansiWeb-method
#' @rdname AnansiWeb
#' @usage NULL
#'
S7::method(show, AnansiWeb) <- function(object) {
    cat(
        class(object),
        " S7 object with ",
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
    cat("Use $ to access: tableX, tableY, dictionary, metadata.")
    invisible(NULL)
}

#' @name dimnames.AnansiWeb
#' @rdname AnansiWeb
#' @method dimnames AnansiWeb
#'
S7::method(dimnames, AnansiWeb) <- function(x) dimnames(x@dictionary)

#' @name dim.AnansiWeb
#' @rdname AnansiWeb
#' @method dim AnansiWeb
#'
S7::method(dim, AnansiWeb) <- function(x) dim(x@dictionary)

#' @name names.AnansiWeb
#' @rdname AnansiWeb
#' @method names AnansiWeb
#'
S7::method(names, AnansiWeb) <- function(x) names(dimnames(x@dictionary))


################################################################################
################################################################################

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
