#' Accessing and modifying information in AnansiWeb S4 class
#' @name AnansiWeb-methods
#' @description `AnansiWeb` supports `$` operator for getting and
#' assigning values.
#'
#' ` dimnames( x ) ` is shorthand for `dimnames( x$dictionary )` and
#' `names( x )` is in turn shorthand for `names( dimnames(x) )`.
#'
#' @returns a specified `AnansiWeb` object.
#' @param ... further arguments.
#'
#' @seealso \itemize{
#' \item [AnansiWeb-class()].
#' \item [weaveWeb()]: for general use.
#' }
#' @importFrom methods slotNames slot slot<-
#' @examples
#' # prepare an AnansiWeb
#' w <- weaveWeb(cpd ~ ko, link = kegg_link())
#'
#' w$dictionary
#'
#' names(w)
#'
NULL

#' @noRd
#' @export
#' @importFrom utils .DollarNames
.DollarNames.AnansiWeb <- function(x, pattern = "") {
    grep(pattern, slotNames(x), value = TRUE)
}

#' @exportMethod $
#' @inheritParams base::`$`
#' @rdname AnansiWeb-methods
#'
setMethod("$", "AnansiWeb", definition = function(x, name) slot(x, name))

#' @exportMethod $<-
#' @importFrom S4Vectors metadata
#' @inheritParams base::`$<-`
#' @rdname AnansiWeb-methods
#'
setReplaceMethod("$", "AnansiWeb", def = function(x, name, value) {
    slot(x, name) <- value
    return(x)
})

#' @export
#' @importClassesFrom S4Vectors Annotated
#' @importFrom S4Vectors metadata
#' @inheritParams S4Vectors::metadata
#' @param simplify `boolean`. If `TRUE` (Default), handles single data.frame
#'     arguments while ensuring compatibility with `S4Vectors` method.
#' @rdname AnansiWeb-methods
#' @importFrom methods slot
#'
setMethod("metadata",
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
#' @rdname AnansiWeb-methods
#'
setReplaceMethod("metadata", "AnansiWeb", def = function(
    x, ..., simplify = TRUE, value) {
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
})

#' @rdname AnansiWeb-methods
#' @name tableY
#' @param x `AnansiWeb`
#' @param ... additional arguments (currently not used).
#' @aliases tableY
#' @export
#'
setMethod("tableY", "AnansiWeb", def = function(x, ...) x@tableY )

#' @rdname AnansiWeb-methods
#' @inheritParams tableY
#' @export
#' @aliases tableX
#'
setMethod("tableX", "AnansiWeb", def = function(x, ...) x@tableX )

#' @rdname AnansiWeb-methods
#' @inheritParams tableY
#' @export
#'
setMethod("dictionary", "AnansiWeb", def = function(x, ...) x@dictionary)

#' @rdname AnansiWeb-methods
#' @name `tableY<-`
#' @inheritParams tableY
#' @param value replacement `matrix` with same number of rows target.
#' @importFrom methods slot<-
#'
setReplaceMethod("tableY", "AnansiWeb", def = function(x, ..., value) {
    x@tableY <- value
    validObject(x)
    x
})

#' @rdname AnansiWeb-methods
#' @export
#' @inheritParams `tableY<-`
#' @aliases `tableX<-`
#'
setReplaceMethod("tableX", "AnansiWeb", def = function(x, ..., value) {
    x@tableX <- value
    validObject(x)
    x
})

#' @rdname AnansiWeb-methods
#' @inheritParams `tableY<-`
#' @export
#'
setReplaceMethod("dictionary", "AnansiWeb", def = function(x, ..., value) {
    x@dictionary <- value
    validObject(x)
    x
})

#' @description `show`: Display the object
#' @importFrom methods show
#' @inheritParams methods::show
#' @rdname AnansiWeb-methods
#' @export
#'
setMethod("show", "AnansiWeb", def = function(object) {
    cat(class(object), " object with ", NROW(object$tableX), " observations:\n",
        "    Tables: ", names(object)[1], " (", NROW(object), " features) and ",
        names(object)[2], " (", NCOL(object), " features)\n",
        sep = ""
    )
    cat("Access content with $ operator. ",
        "Collapse with as.list().",
        sep = ""
    )
    invisible(NULL)
})

#' @rdname AnansiWeb-methods
#' @inheritParams base::dimnames
#' @export
#'
setMethod(
    "dimnames", "AnansiWeb",
    function(x) dimnames(x@dictionary)
)

#' @rdname AnansiWeb-methods
#' @inheritParams base::dim
#' @export
#'
setMethod(
    "dim", "AnansiWeb",
    function(x) dim(x@dictionary)
)

#' @rdname AnansiWeb-methods
#' @inheritParams base::names
#' @export
#'
setMethod("names", "AnansiWeb", function(x) names(dimnames(x@dictionary)))

#' @noRd
#'
tell_F <- function(tale) {
    if (is(tale, "anansiTale")) {
        return(tale@f.values)
    }
    if (is.list(tale)) {
        return(lapply(tale, tell_F))
    }
}

#' @noRd
#'
tell_T <- function(tale) {
    if (is(tale, "anansiTale")) {
        return(tale@t.values)
    }
    if (is.list(tale)) {
        return(lapply(tale, tell_T))
    }
}

#' @noRd
#'
tell_P <- function(tale) {
    if (is(tale, "anansiTale")) {
        return(tale@p.values)
    }
    if (is.list(tale)) {
        return(lapply(tale, tell_P))
    }
}


#' @noRd
#'
tell_e <- function(tale) {
    if (is(tale, "anansiTale")) {
        return(tale@estimates)
    }
    if (is.list(tale)) {
        return(lapply(tale, tell_e))
    }
}


#' @noRd
#'
tell_df1 <- function(tale) {
    if (is(tale, "anansiTale")) {
        return(tale@df[1])
    }
    if (is.list(tale)) {
        return(lapply(tale, tell_df1))
    }
}

#' @noRd
#'
tell_df2 <- function(tale) {
    if (is(tale, "anansiTale")) {
        return(tale@df[2])
    }
    if (is.list(tale)) {
        return(lapply(tale, tell_df2))
    }
}

#' @noRd
#'
tell_dfr <- function(tale) {
    if (is(tale, "anansiTale")) {
        return(tale@df[3])
    }
    if (is.list(tale)) {
        return(lapply(tale, tell_dfr))
    }
}

################################################################################
################################################################################

#' Is this a data.frame with exactly two columns that are named?
#' @noRd
validWeb <- function(x) {
    y_names <- identical(rownames(x), colnames(x$tableY))
    x_names <- identical(colnames(x), colnames(x$tableX))
    s_names <- identical(rownames(x$tableY), rownames(x$tableX))
    meta_dim <- any(
        NROW(x$metadata) == NROW(x$tableY),
        prod(dim(x$metadata)) <= 1
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
