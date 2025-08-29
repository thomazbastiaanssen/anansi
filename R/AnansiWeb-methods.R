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
#' @param x input, `AnansiWeb`.
#' @param ... further arguments for methods
#' @usage
#' ## Accessors
#' dimnames(x)
#' dim(x)
#' names(x)
#'
#' ## Coercion
#' asMAE(x)
#' as.list(x, ...)
#' as.data.frame(x, row.names = NULL, optional = FALSE, ...)
#'
#' ## Utilities on feature pairs
#' pairs(x, ...)
#' getFeaturePairs(
#'     x, which = NULL, with.metadata = FALSE, ...
#' )
#' pairwiseApply(
#'     x,
#'     FUN,
#'     MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE
#' )
#'
#'
#' @seealso \itemize{
#' \item [weaveWeb()]: for general use.
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
#'
#' # Extract data.frames in pairs (only show first)
#' getFeaturePairs(web)[1L]
#'
#' pairwiseApply(
#'     FUN = function(x, y) cor(x, y),
#'     web
#' )
NULL

#' @name AnansiWeb
#' @importFrom methods show
#' @importMethodsFrom methods show
#' @aliases show,anansi::AnansiWeb-method
#' @rdname AnansiWeb
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

#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases dimnames.anansi::AnansiWeb
#' @export
#' @usage NULL
#'
S7::method(dimnames, AnansiWeb) <- function(x) dimnames(x@dictionary)

#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases dim.anansi::AnansiWeb
#' @export
#' @usage dim(x)
#'
S7::method(dim, AnansiWeb) <- function(x) dim(x@dictionary)

#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases names.anansi::AnansiWeb
#' @export
#' @usage names(x)
#'
S7::method(names, AnansiWeb) <- function(x) names(dimnames(x@dictionary))

#' @name AnansiWeb
#' @rdname AnansiWeb
#' @importFrom graphics pairs
#' @aliases pairs pairs,anansi::AnansiWeb-method pairs.anansi::AnansiWeb
#' @returns  a two-column array index, corresponding to i,j coordinates in
#'   matrix notation.
#' @export
#' @usage NULL
#'
S7::method(pairs, AnansiWeb) <- function(x, ...) Matrix::which(
    x, arr.ind = TRUE, useNames = FALSE
    )

#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases pairwiseApply pairwiseApply.anansi::AnansiWeb
#' @param FUN a function with at least two arguments. The variables `x` and `y`,
#'     in order, refer to the corresponding values of feature pairs in `tableX`
#'     and `tableY`.
#' @param MoreArgs,SIMPLIFY,USE.NAMES see ?base::mapply
#' @export
#' @usage NULL
#'
S7::method(pairwiseApply, AnansiWeb) <- function(
        X, FUN, MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE
        ) {

        tY <- as.data.frame.matrix(X@tableY, make.names = FALSE)
        tX <- as.data.frame.matrix(X@tableX, make.names = FALSE)
        wh <- Matrix::which(X@dictionary, arr.ind = TRUE, useNames = FALSE)

        out <- base::.mapply(
            FUN,
            dots = list(
                x = tX[wh[, 2L]],
                y = tY[wh[, 1L]]
            ),
            MoreArgs
        )
        if (USE.NAMES) {
            names(out) <- paste0(colnames(tX)[wh[, 2L]], colnames(tY)[wh[, 1L]])
        }

        if (SIMPLIFY) {
            out <- simplify2array(out)
        }
        return(out)
    }

#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases getFeaturePairs getFeaturePairs,anansi::AnansiWeb-method
#' @importFrom Matrix which
#' @param which `integer matrix`, indicating pair positions in `x@tableY` and
#'     `x@tableX`, respectively. If `NULL` (default):
#'     `Matrix::which(x@dictionary, TRUE)`.
#' @param with.metadata `Logical scalar` whether to append metadata to output
#' @return A list of data.frames with the paired data
#' @usage NULL
#' @export
#'
S7::method(getFeaturePairs, AnansiWeb) <-  function(x, ..., which = NULL, with.metadata = FALSE) {
    if (is.null(which)) {
        which <- Matrix::which(x@dictionary, arr.ind = TRUE, useNames = FALSE)
    }
    tX <- x@tableX
    tY <- x@tableY
    xnames <- colnames(tX)
    ynames <- colnames(tY)
    if (!with.metadata) {
        return(
            lapply(seq_len(NROW(which)), FUN = function(z) {
                cbind(
                    tY[, which[z, 1L], drop = FALSE],
                    tX[, which[z, 2L], drop = FALSE]
                )
            })
        )
    } else {
        metadata <- x@metadata
        return(
            lapply(seq_len(NROW(which)), FUN = function(z) {
                cbind(
                    tY[, which[z, 1L], drop = FALSE],
                    tX[, which[z, 2L], drop = FALSE],
                    metadata
                )
            })
        )
    }
}

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
