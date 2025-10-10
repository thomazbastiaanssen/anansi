#' Methods for pairwise operations on AnansiWeb objects
#' @description
#' Methods to run pairwise analysis on multi-modal data.
#'
#' @name AnansiWeb-pairwise
#' @rdname AnansiWeb-pairwise
#' @param X,x input, AnansiWeb object
#' @param ... additional arguments
#' @param FUN a function with at least two arguments. The variables `x` and `y`,
#'     in order, refer to the corresponding values of feature pairs in `tableX`
#'     and `tableY`.
#' @param MoreArgs,SIMPLIFY,USE.NAMES see ?base::mapply
#' @param which `integer matrix`, indicating pair positions in `x@tableY` and
#'     `x@tableX`, respectively. If `NULL` (default):
#'     `Matrix::which(x@dictionary, TRUE)`.
#' @param with.metadata `Logical scalar` whether to append metadata to output
#' @aliases pairs.AnansiWeb pairs.anansi::AnansiWeb pairwiseApply.AnansiWeb
#' @aliases pairwiseApply.anansi::AnansiWeb
#' @returns
#' pairs: A two-column array index, corresponding to i,j coordinates in matrix
#' notation. \cr
#' getFeaturePairs: A list of data.frames with the paired data \cr
#'
#'

#' @examples
#' web <- randomWeb()
#' # Extract data.frames in pairs (only show first)
#' getFeaturePairs(web)[1L]
#'
#' pairs(x = web)
#'
#' getFeaturePairs(x = web, which = NULL, with.metadata = FALSE)
#'
#' pairwiseApply(
#'     X = web,
#'     FUN = function(x, y) cor(x, y),
#'     MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE
#' )
#'
NULL

#' @export
#' @importFrom graphics pairs
S7::method(pairs, AnansiWeb) <- function(x, ...) `pairs.anansi::AnansiWeb`(x)

#' @export
`pairs.anansi::AnansiWeb` <- function(x, ...) {
    Matrix::which(
        x@dictionary,
        arr.ind = TRUE, useNames = FALSE
    )
}

#' @importFrom Matrix which
#'
S7::method(getFeaturePairs, AnansiWeb) <- function(
        x, which = NULL, with.metadata = FALSE
) {
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


#' @export
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
