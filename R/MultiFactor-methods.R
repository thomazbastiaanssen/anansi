#' S4 Methods for MultiFactor
#' @name MultiFactor-methods
#' @description
#' a set of methods to work with `MultiFactor`, a special type of list.
#'
NULL


#' S4 Methods for MultiFactor
#' @description
#' `getEdgeList`: Return a data frame in edge list format.
#' @rdname MultiFactor-methods
#' @export
#'
setMethod(
    "getEdgeList", "MultiFactor",
    function(x) as.data.frame(do.call(rbind, names(x)))
)


#' S4 Methods for MultiFactor
#' @description `dim`: Display a vector of dims (n objects, n ids).
#' @rdname MultiFactor-methods
#' @export
#'
setMethod("dim", "MultiFactor", function(x) {
    dim(x@map)
})

#' S4 Methods for MultiFactor
#' @description `names`: Display a vector of dims (n objects, n ids).
#' @rdname MultiFactor-methods
#' @export
#'
setMethod("names", "MultiFactor", function(x) {
    `names<-`(lapply(x@index, names), rownames(x))
})

#' S4 Methods for MultiFactor
#' @description `dimnames`: Display a vector of dims (n objects, n ids).
#' @rdname MultiFactor-methods
#' @export
#'
setMethod("dimnames", "MultiFactor", function(x) {
    dimnames(x@map)
})

#' S4 Methods for MultiFactor
#' @description `[`: Subset based on `[rownames(),colnames(x)]` indices.
#' @param drop Whether to return a `list` (Default) or `MultiFactor`.
#' @export
#'
setMethod("[", c("MultiFactor", "ANY", "ANY"), definition = function(
    x, i, j, ..., drop = TRUE) {
    if (missing(i) && missing(j)) {
        return(x)
    }
    d <- dictionary(x)
    l <- levels(x)
    x <- x@index
    if (!missing(i)) ii <- rownames(d[i, , drop = FALSE])
    if (!missing(j)) jj <- colnames(d[, j, drop = FALSE])

    if (missing(i)) {
        ii <- rowsWithCol(d, jj, FALSE)
        x <- lapply(x[ii], `[`, i = jj)
    } else if (missing(j)) {
        x <- x[ii]
    } else {
        x <- lapply(x[ii], `[`, i = jj)
    }

    if (drop) {
        return(x)
    }

    MultiFactor(x, levels = l)
})

#' S4 Methods for MultiFactor
#' @description `[<-`: Replace based on `[rownames(),colnames(x)]` indices
#' @export
#'
setReplaceMethod("[", c("MultiFactor", "ANY", "ANY", "list"), def = function(
    x, i, j, ..., value) {
    if (missing(i) && missing(j)) {
        return(value)
    }
    d <- dictionary(x)
    if (!missing(i)) ii <- rownames(d[i, , drop = FALSE])
    if (!missing(j)) jj <- colnames(d[, j, drop = FALSE])

    if (missing(j)) {
        x@index[ii] <- value
        validObject(x)
        return(x)
    }

    if (missing(i)) {
        ii <- rowsWithCol(d, jj, names = TRUE)
    }

    for (i in ii) {
        for (j in jj) {
            x@index[[i]][, j] <- value[[i]][, j]
        }
    }
    validObject(x)
    (x)
})



#' @export
#' @description
#' Extract all elements that contain a given feature type.
#' @details
#' if argument `i` is a vector, return all elements that contain at least that
#' exact combination.
#' @rdname MultiFactor-methods
#' @inheritParams base::Extract
#'
setMethod("[[", c("MultiFactor", "ANY"), function(x, i, ...) {
    d <- x@map
    # If i can't index d, return NULL
    if (!all(i %in% colnames(d))) {
        if (anyNA(colnames(d)[i]) || length(colnames(d)[i]) != length(i)) {
            return(NULL)
        }
    }
    # Otherwise, return selected elements.
    ii <- rowsWithCol(d, i, FALSE)
    x[ii]
})

#' @export
#' @rdname MultiFactor-methods
#' @description `[[<-` should be used to replace, not add data.
#'
setReplaceMethod("[[", c("MultiFactor", "ANY", "ANY"), function(
    x, i, ..., value) {
    d <- x@map
    # If i can't index d, stop. Appending not supported through `[[<-`.
    if (!all(i %in% colnames(d))) {
        if (anyNA(colnames(d)[i]) || length(colnames(d)[i]) != length(i)) {
            stop("No levels corresponding to `i` found in MultiFactor. ")
        }
    }
    ii <- rowsWithCol(d, i, FALSE)
    x@index[ii] <- value
    validObject(x)
    x
})

ff <- function(x, ...){
args <- list(...)
if (!length(args)) { return(x) }

index <- c(x, ...)
index
}
#' @export
#' @rdname MultiFactor-methods
#' @description Combine two or more `MultiFactor` objects.
#'
setMethod("c", "MultiFactor", function(x, ...) {

    args <- list(...)
    if (!length(args)) { return(x) }

    ind <- Reduce(mergeMultiFactorInds, list(x@index, ...))

    ind <- checkMergers(ind, TRUE)
    lvs <- Reduce(mergeMultiFactorLvs, list(levels(x), ...))
    MultiFactor(ind, levels = lvs)

    return(x)
})

#' @noRd
#' @param `MultiFactor@index` from first `MultiFactor` in `c()` Method.
#' @param y a second `MultiFactor`
#' @returns Merged index.
mergeMultiFactorInds <- function(x, y) c(x, y@index)

#' @noRd
#' @param x `levels(MultiFactor)` from first `MultiFactor` in `c()` Method.
#' @param y a second `MultiFactor`.
#' @returns Merged levels
mergeMultiFactorLvs <- function(x, y) {
    y <- levels(y)
    i <- intersect(names(x), names(y))
    x[i] <- union(x[i], y[i])
    return(c(x, y[! names(y) %in% i]))
    }

#' S4 Methods for MultiFactor
#' @description `show`: Display the object
#' @importFrom methods show
#' @importFrom Matrix sparseMatrix printSpMatrix
#' @inheritParams methods::show
#' @rdname MultiFactor-methods
#' @export
#'
setMethod("show", "MultiFactor", function(object) {
    cat("A list of class ", class(object), ",\n    ",
        NCOL(object), " feature types across ", NROW(object),
        " edge lists.\n\n",
        sep = ""
    )
    printSpMatrix(object@map)

    cat("\nValues represent unique feature names in that edge list.\n\n",
        "Levels:\n\n",
        sep = ""
    )
    id_w <- max(nchar(colnames(object)))
    nm_w <- max(nchar(nlevels(object)))
    for (id in colnames(object)) {
        num_lvs <- length(levels(object)[[id]])
        cat(format(id, width = id_w), " : ",
            format(num_lvs, width = nm_w), " Levels: ",
            sep = ""
        )

        if (num_lvs > 4L) {
            cat(
                levels(object)[[id]][1], levels(object)[[id]][2], "...",
                levels(object)[[id]][num_lvs], "\n",
                sep = " "
            )
        } else {
            cat(levels(object)[[id]], "\n", sep = " ")
        }
    }
    invisible(NULL)
})

#' @rdname MultiFactor-methods
#' @description Analogous to `factors`. `unfactor(MultiFactor)` returns a named
#'     list with character data frames with same dimensions as input.
#' @importMethodsFrom S4Vectors unfactor
#' @inheritParams S4Vectors::unfactor
#' @returns A named character list
#' @export
#'
setMethod("unfactor", "MultiFactor", function(x) {
    lv <- levels(x)
    ns <- rownames(x)
    x <- x@index

    x[] <- lapply(x, function(df) {
        for (id in names(df)) {
            df[, id] <- lv[[id]][df[, id]]
        }
        return(df)
    })

    return(x)
})

#' @rdname MultiFactor-methods
#' @description Analogous to `factors`. `droplevels(MultiFactor)` returns a
#'     `MultiFactor` with unused levels removed.
#' @importMethodsFrom S4Vectors droplevels
#' @inheritParams base::droplevels
#' @param exclude `NULL` or `Named character list` of similar structure as
#'     `levels(MultiFactor)`. Which levels to drop from output.
#' @param select `NULL` or `Named character list` of similar structure as
#'     `levels(MultiFactor)`. Which levels to keep in output.
#' @details Only one of `select` and `exclude` should be provided, as they are
#'     each others complement.
#' @returns A MultiFactor
#' @export
#'
droplevels.MultiFactor <- function(x, exclude = NULL, select = NULL, ...) {
    stopifnot(
        "Only one of 'exclude' and 'select' may be provided" =
            sum(is.null(exclude), is.null(select)) > 0L
    )
    stopifnot("'x' is not a MultiFactor." = is(x, "MultiFactor"))
    # Section 1. Trimming the indices by user input
    lvs <- levels(x)
    d <- dictionary(x)
    if (!is.null(exclude)) {
        stopifnot(
            "`'exclude' must be a named list of character vectors ." =
                is.list(exclude) && any(names(exclude) %in% names(lvs))
        )
        # Names not mentioned will be left alone.
        jj <- intersect(names(exclude), names(lvs))
        for (j in jj) {
            ex.ind <- match(exclude[[j]], lvs[[j]], nomatch = 0L)
            ii <- rowsWithCol(d, j)
            for (i in ii) {
                x@index[[i]] <- x@index[[i]][!x@index[[i]][, j] %in% ex.ind, ]
            }
        }
    } else if (!is.null(select)) {
        stopifnot(
            "'select' argument must be a named list of character vectors." =
                is.list(select) && any(names(select) %in% names(lvs))
        )
        jj <- intersect(names(select), names(lvs))
        for (j in jj) {
            ex.ind <- match(select[[j]], lvs[[j]], nomatch = 0L)
            ii <- rowsWithCol(d, j)
            for (i in ii) {
                x@index[[i]] <- x@index[[i]][x@index[[i]][, j] %in% ex.ind, ]
            }
        }
    }
    # Section 2. Trimming the levels by indices. .
    for (lv in names(lvs)) {
        # Loop over all cols. First determine which rows are relevant per col/type
        rs <- rowsWithCol(d, lv, names = TRUE)
        x_index <- lapply(x@index[rs], `[[`, lv)
        # Get unique feature names in that type and are within levels.
        x_tot <- unique(unlist(x_index, use.names = FALSE))
        keep_ix <- which(seq_along(lvs[[lv]]) %in% x_tot)
        # Keep levels that show up in data
        x@levels[[lv]] <- lvs[[lv]][keep_ix]
        # Update indices to reflect fewer level names.
        for (r in rs) {
            x@index[[r]][, lv] <- match(x_index[[r]], table = keep_ix)
        }
    }
    x@map <- mapMultiFactor(x@index, mode = "counts")
    validObject(x)
    return(x)
}

#' @rdname MultiFactor-methods
#' @export
#'
setMethod("droplevels", "MultiFactor", function(x, ...) {
    droplevels.MultiFactor(x, ...)
})


# mergeROWS()

#' S3/S4 combo for levels.
#' @export
#' @description
#' get object levels
#' @returns a named list of character vectors.
#' @rdname MultiFactor-methods
#'
levels.MultiFactor <- function(x) x@levels

#' @export
#' @rdname MultiFactor-methods
setMethod("levels", "MultiFactor", levels.MultiFactor)

#' @export
#' @rdname MultiFactor-methods
#' @param value a replacement character vector of suitable dimensions.
#'
setReplaceMethod(
    "levels", "MultiFactor",
    function(x, value) {
        x@levels <- value
        x
    }
)

#' @export
#' @description
#' get object map
#' @param x `MultiFactor` object
#' @returns a named sparse biadjacency matrix of dimensions (`dimnames(x)`)
#' @rdname MultiFactor-methods
#'
setMethod("dictionary", "MultiFactor", function(x, ...) x@map)

#' @export
#' @rdname MultiFactor-methods
setReplaceMethod(
    "dictionary", "MultiFactor",
    function(x, ..., value) {
        x@map <- value
        validObject(x)
        x
    }
)


#' S4 Methods for MultiFactor
#' @rdname MultiFactor-methods
#' @details
#' `subset`: For `MultiFactor objects`, sub-setting is only applied
#' to data frames compatible with the expression. The rest are returned
#' unaltered. Modeled after `subset()`.
#'
#' @param subset
#' `logical expression` indicating rows to keep. Must contain variables
#' found as column names.
#' @param select `expression`. Which column names to consider. If missing
#' (Default), consider all column names.
#' @inheritParams BiocGenerics::subset
#' @importMethodsFrom BiocGenerics subset
#' @export
#' @seealso [BiocGenerics::subset()].
#' [weaveWeb()] for the AnansiWeb constructor functions that
#' take link data frames.
#' @examples
#' # prep input
#' l <- asMultiFactor(kegg_link())
#'
#' # Sub-setting is only performed on data frames that contain the arguments
#' str(subset(x = l, cpd %in% c("C00001", "C00002")))
#'
#' # Several data frames at the same time:
#' subset(x = l, ec %in% c("1.2.3.4", "4.3.2.1"))
#'
setMethod("subset", "MultiFactor", function(x, subset, select, ...) {
    validObject(x)
    x.names <- names(x)
    # PART I: SUBSETTING
    if (!missing(subset)) {
        subset <- substitute(subset)
        sub.vars <- all.vars(subset)
        # Select those data frames where all terms are mentioned
        sub.ind <- unlist(lapply(x.names, function(y) all(sub.vars %in% y)))
        # Subset them
        x[sub.ind] <- lapply(x[sub.ind], function(y) {
            r <- eval(subset, y, parent.frame())
            return(y[r, ])
        })
    }
    # Return now if only one df.
    if (length(x) == 1L) {
        return(x)
    }

    # PART II: SELECTING
    if (missing(select)) {
        id.vec <- unlist(x.names, use.names = FALSE)
        id.share <- id.vec[duplicated(id.vec)]
        sel.vars <- id.share
    } else {
        sel.vars <- all.vars(substitute(select))
    }
    for (v in sel.vars) {
        # Select those data frames where all terms are mentioned
        s.ind <- unlist(lapply(x.names, function(y) v %in% y))
        sel.obj <- x[s.ind]
        keep <- Reduce(intersect, lapply(sel.obj, function(df) df[, v]))
        # Filter feature ids in each df to only include universally shared ones.
        x[s.ind] <- lapply(sel.obj, function(df) {
            return(df[df[[v]] %in% keep, ])
        })
    }
    return(x)
})


#' @param d `MultiFactor@map`
#' @param id `Character or Integer scalar`. Selects column(s) of `d`.
#' @param names Whether to return characters (Default) or integer indices.
#' @returns A vector indicating which elements of `MultiFactor` contain `id`.
#' @importFrom Matrix rowSums
#' @noRd
#' @description Helper function for `MultiFactor` to get names or indices of
#' data frames that contain an id column
#'
rowsWithCol <- function(d, id, names = TRUE) {
    rowInds <- which(Matrix::rowSums(d[, id, drop = FALSE] > 0L) == length(id))
    if (length(rowInds) == 0L) {
        return(NULL)
    }
    if (names) {
        rowInds <- rownames(d)[rowInds]
    }
    return(rowInds)
}

#' @noRd
#' @description `rowsWithCol` but returns union rather than intersect.
rowsInCol <- function(d, id, names = TRUE) {
    rowInds <- which(Matrix::rowSums(d[, id, drop = FALSE] > 0L) > 0L)
    if (length(rowInds) == 0L) {
        return(NULL)
    }
    if (names) {
        rowInds <- rownames(d)[rowInds]
    }
    return(rowInds)
}

#' @param d `MultiFactor@map`
#' @param id `Character or Integer vector`. Selects row(s) of `d`.
#' @param names Whether to return characters (Default) or integer indices.
#' @returns A vector indicating which feature types are in element `id`.
#' @importFrom Matrix colSums
#' @noRd
#' @description Helper function for `MultiFactor` to get names or indices of
#'     features contained in a given data frame element of `MultiFactor`.
#'
colsWithRow <- function(d, id, names = TRUE) {
    colInds <- which(Matrix::colSums(d[id, , drop = FALSE] > 0L) == length(id))
    if (length(colInds) == 0L) {
        return(NULL)
    }
    if (names) {
        colInds <- colnames(d)[colInds]
    }
    return(colInds)
}

#' @noRd
#' @description `colsWithRow` but returns union rather than intersect.
colsInRow <- function(d, id, names = TRUE) {
    colInds <- which(Matrix::colSums(d[id, , drop = FALSE] > 0L) > 0L)
    if (length(colInds) == 0L) {
        return(NULL)
    }
    if (names) {
        colInds <- colnames(d)[colInds]
    }
    return(colInds)
}
