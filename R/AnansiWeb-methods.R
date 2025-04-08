#' AnansiWeb S4 container class
#' @name AnansiWeb
#' @description
#' `AnansiWeb` is an S4 class containing two feature tables as well as a
#' dictionary to link them. `AnansiWeb` is the main container that will
#' hold your input data throughout the `anansi` pipeline.
#'
#' Typical use of the `anansi` package will involve generating an `AnansiWeb`
#' object using the `weaveWeb()` function.
#'
#' The function `AnansiWeb()` constructs an `AnansiWeb` object from two
#' feature tables and an adjacency matrix.
#'
#' @usage
#' ## Accessors
#' \S4method{dimnames}{AnansiWeb}(x)
#' \S4method{dim}{AnansiWeb}(x)
#' \S4method{names}{AnansiWeb}(x)
#'
#' \S4method{tableY}{AnansiWeb}(x, ...)
#' \S4method{tableY}{AnansiWeb}(x, ...) <- value
#' \S4method{tableX}{AnansiWeb}(x, ...)
#' \S4method{tableX}{AnansiWeb}(x, ...) <- value
#' \S4method{dictionary}{AnansiWeb}(x, ...)
#' \S4method{dictionary}{AnansiWeb}(x, ...) <- value
#' \S4method{metadata}{AnansiWeb}(x, simplify = TRUE, ...)
#' \S4method{metadata}{AnansiWeb}(x, simplify = TRUE, ...) <- value
#'
#' ## Coercion
#' asMAE(x)
#' \S4method{as.list}{AnansiWeb}(x, ...)
#' \S4method{as.data.frame}{AnansiWeb}(
#'     x, row.names = NULL, optional = FALSE, ...
#'     )
#'
#' ## Utilities on feature pairs
#' \S4method{which}{AnansiWeb}(x, arr.ind = TRUE, useNames = FALSE)
#' \S4method{getFeaturePairs}{AnansiWeb}(
#'     x, which = NULL, with.metadata = FALSE, ...
#' )
#' \S4method{mapply}{AnansiWeb}(
#'     FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE
#'     )
#'
#' @param x,object an `AnansiWeb` object on which a method will be applied.
#'
#' @seealso \itemize{
#' \item [weaveWeb()]: for general use.
#' }
#' @importFrom methods slotNames slot slot<-
#' @aliases tableY `tableY<-`
#' @examples
#'
#' # Methods for AnansiWeb
#' dimnames(web)
#' dim(web)
#' names(web)
#'
#' tableX(web)
#' tableY(web)
#' dictionary(web)
#'
#' # Assign some random metadata
#' metadata(web) <- data.frame(
#'     id = row.names(tableY(web)),
#'     a = rnorm(36),
#'     b = sample(c("a", "b"), 36, TRUE),
#'     row.names = "id"
#' )
#' metadata(web)
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
#' mapply(
#'     FUN = function(x, y) cor(x, y),
#'     web
#' )
#'
NULL

#' @export
#' @importClassesFrom S4Vectors Annotated
#' @importMethodsFrom S4Vectors metadata
#' @param simplify `boolean`. If `TRUE` (Default), handles single data.frame
#'     arguments while ensuring compatibility with `S4Vectors` method.
#' @rdname AnansiWeb
#' @aliases dictionary metadata,AnansiWeb-method
#' @importFrom methods slot
#' @usage NULL
#'
setMethod(
    "metadata",
    signature = c(x = "AnansiWeb"),
    definition = function(x, simplify = TRUE, ...) {
        m <- x@metadata
        if (simplify && "metadata" %in% names(m)) {
            return(m[["metadata"]])
        }
        return(m)
    }
)

#' @export
#' @importMethodsFrom S4Vectors "metadata<-"
#' @importFrom methods slot<-
#' @aliases metadata<-,AnansiWeb-method
#' @rdname AnansiWeb
#' @usage NULL
#'
setReplaceMethod(
    "metadata",
    "AnansiWeb",
    def = function(
        x,
        ...,
        simplify = TRUE,
        value
    ) {
        if (simplify && inherits(value, "data.frame")) {
            x@metadata[["metadata"]] <- as.data.frame(value)
            return(x)
        }
        if (!is.list(value)) {
            stop("replacement 'metadata' value must be a list")
        }
        if (!length(value)) {
            names(value) <- NULL
        } # instead of character()
        x@metadata <- value
        validObject(x)
        x
    }
)

#' @rdname AnansiWeb
#' @param ... additional arguments (currently not used).
#' @aliases tableY tableY,AnansiWeb-method
#' @export
#' @usage NULL
#'
setMethod("tableY", "AnansiWeb", def = function(x, ...) x@tableY)

#' @rdname AnansiWeb
#' @export
#' @aliases tableX tableX,AnansiWeb-method
#' @usage NULL
#'
setMethod("tableX", "AnansiWeb", def = function(x, ...) x@tableX)

#' @rdname AnansiWeb
#' @aliases `dictionary` dictionary,AnansiWeb-method
#' @export
#' @usage NULL
#'
setMethod("dictionary", "AnansiWeb", def = function(x, ...) x@dictionary)

#' @rdname AnansiWeb
#' @aliases tableY<- tableY<-,AnansiWeb-method
#' @param value replacement `matrix` with same number of rows target.
#' @usage NULL
#'
setReplaceMethod("tableY", "AnansiWeb", def = function(x, ..., value) {
    x@tableY <- value
    validObject(x)
    x
})

#' @rdname AnansiWeb
#' @export
#' @aliases tableX<- tableX<-,AnansiWeb-method
#' @usage NULL
#'
setReplaceMethod("tableX", "AnansiWeb", def = function(x, ..., value) {
    x@tableX <- value
    validObject(x)
    x
})

#' @rdname AnansiWeb
#' @aliases dictionary<- dictionary<-,AnansiWeb-method
#' @export
#' @usage NULL
#'
setReplaceMethod("dictionary", "AnansiWeb", def = function(x, ..., value) {
    x@dictionary <- value
    validObject(x)
    x
})

#' @importFrom methods show
#' @rdname AnansiWeb
#' @export
#'
setMethod("show", "AnansiWeb", def = function(object) {
    cat(
        class(object),
        " S4 object with ",
        NROW(tableX(object)),
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
    cat("Accessors: tableX(), tableY(), dictionary(), metadata().")
    invisible(NULL)
})

#' @rdname AnansiWeb
#' @export
#' @usage NULL
#'
setMethod(
    "dimnames",
    "AnansiWeb",
    function(x) dimnames(x@dictionary)
)

#' @rdname AnansiWeb
#' @export
#' @usage NULL
#'
setMethod(
    "dim",
    "AnansiWeb",
    function(x) dim(x@dictionary)
)

#' @rdname AnansiWeb
#' @export
#' @usage NULL
#'
setMethod("names", "AnansiWeb", function(x) names(dimnames(x@dictionary)))

#' @rdname AnansiWeb
#' @aliases which,AnansiWeb-method
#' @importMethodsFrom BiocGenerics which
#' @param arr.ind,useNames See ?base::which. `AnansiWeb` default returns a
#'     two-column array index.
#' @export
#' @usage NULL
#'
setMethod(
    "which",
    signature = c(x = "AnansiWeb"),
    function(x, arr.ind = TRUE, useNames = FALSE) {
        which.AnansiWeb(x, arr.ind, useNames)
    }
)

#' @noRd
#' @importMethodsFrom Matrix which
#'
which.AnansiWeb <- function(x, arr.ind = TRUE, useNames = FALSE) {
    Matrix::which(x@dictionary, arr.ind, useNames)
}

#' @rdname AnansiWeb
#' @aliases mapply,AnansiWeb-method
#' @importMethodsFrom BiocGenerics mapply
#' @param FUN a function with at least two arguments. The variables `x` and `y`,
#'     in order, refer to the corresponding values of feature pairs in `tableX`
#'     and `tableY`.
#' @param MoreArgs,SIMPLIFY,USE.NAMES see ?base::mapply
#' @export
#' @usage NULL
#'
setMethod(
    "mapply",
    signature = c(... = "AnansiWeb"),
    function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE) {
        tY <- as.data.frame.matrix(tableY(...), make.names = FALSE)
        tX <- as.data.frame.matrix(tableX(...), make.names = FALSE)
        wh <- which.AnansiWeb(...)

        out <- .mapply(
            FUN,
            dots = list(x = tX[wh[, 2L]], y = tY[wh[, 1L]]),
            MoreArgs
        )
        if (SIMPLIFY) {
            out <- simplify2array(out)
        }
        return(out)
    }
)

#' @rdname AnansiWeb
#' @aliases getFeaturePairs getFeaturePairs,AnansiWeb-method
#' @importFrom Matrix which
#' @param which `integer matrix`, indicating pair positions in `tableY(x)` and
#'     `tableX(x)`, respectively. If `NULL` (default):
#'     `Matrix::which(dictionary(x), TRUE)`.
#' @param with.metadata `Logical scalar` whether to append metadata to output
#' @return A list of data.frames with the paired data
#' @usage NULL
#' @export
#'
setMethod(
    getFeaturePairs,
    "AnansiWeb",
    function(x, which = NULL, with.metadata = FALSE, ...) {
        getFeaturePairs.AnansiWeb(x, which, with.metadata)
    }
)

#' @rdname AnansiWeb
#' @noRd
getFeaturePairs.AnansiWeb <- function(x, which = NULL, with.metadata = FALSE) {
    if (is.null(which)) {
        which <- which(x)
    }
    tX <- tableX(x)
    tY <- tableY(x)
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
        metadata <- metadata(x)
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
    y_names <- identical(rownames(x), colnames(tableY(x)))
    x_names <- identical(colnames(x), colnames(tableX(x)))
    s_names <- identical(rownames(tableY(x)), rownames(tableX(x)))
    meta_dim <- any(
        NROW(metadata(x)) == NROW(tableY(x)),
        prod(dim(metadata(x))) <= 1
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
