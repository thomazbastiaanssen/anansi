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
#'# Setup
#' web <- randomWeb(n_samp = 36)
#'
#' # Accessors
#' dimnames(web)
#' dim(web)
#' names(web)
#'
#' # Getters and setters: `@` and `@<-`
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
#' # Coerce to list
#' weblist <- as.list(web)
#'
#' # Coerce to Data.frame
#' webdf <- as.data.frame(web)
#'
#' # Coerce to MultiAssayExperiment
#' mae <- asMAE(web)
NULL

#' @importFrom methods show
#' @importMethodsFrom methods show
#' @export
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

#' @export
#'
S7::method(dimnames, AnansiWeb) <- function(x) dimnames(x@dictionary)

#' @export
#'
S7::method(dim, AnansiWeb) <- function(x) dim(x@dictionary)

#' @export
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
