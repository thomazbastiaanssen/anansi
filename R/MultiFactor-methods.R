#' MultiFactor S7 container class
#' @name MultiFactor
#' @rdname MultiFactor
#' @description
#' `MultiFactor` is an S4 class to organize and manage multiple sets of factors,
#' for instance when tracing or converting feature IDs across databases. Methods
#' for `MultiFactor` aim to follow `factor` behaviour.
#'
#' @details
#' The most straightforward way to construct a `MultiFactor` object is as a
#' named list of named data.frames. The columns of the data.frames indicate the
#' category of factor in that column.
#'
#' A `MultiFactor` object presents itself similar to a `data.frame`, in the
#' sense that level types can be called as columns and individual data.frame
#' components can be called as rows.
#'
#' @usage
#' ## Accessors
#' dimnames(x)
#' dim(x)
#' names(x)
#'
#' ## Indexing
#' x[...]
#' x[...] <- value
#' x[[...]]
#' x[[...]] <- value
#'
#' ## Factor manipulation
#' levels(x)
#' unfactor(x)
#' droplevels(x, exclude = NULL, select = NULL)
#'
#' ## Coercion
#' as.list(x, ..., use.names = TRUE)
#' @examples
#' x <- MultiFactor(kegg_link())
#' x
#' dimnames(x)
#' levels(x)
#'
#' @param x `MultiFactor` on which the method should be applied, or, in
#'     case of the constructor `MultiFactor()`, a named `list` of data.frames
#'     with two named columns each, where elements that share a row indicates
#'     thet are adjacent.
#' @param ... `i,j` indices specifying elements to extract or replace. Indices are
#'     numeric or character vectors or empty (missing) or NULL. Numeric values
#'     are coerced to integer or whole numbers as by as.integer or for large
#'     values by trunc (and hence truncated towards zero). Character vectors
#'     will be matched to the names of the object.
#' @param value Replacement value, typically of same type as that which is to be
#'     replaced.
NULL

#' @name MultiFactor
#' @rdname MultiFactor
#' @aliases getEdgeList
#' @usage getEdgeList(x)
#' @export
#'
S7::method(getEdgeList, MultiFactor) <- function(x) {
    as.data.frame(do.call(rbind, base::names(x)))
    }

#' @name MultiFactor
#' @rdname MultiFactor
#' @export
S7::method(dim, MultiFactor) <- function(x) dim(x@map)



#' @name MultiFactor
#' @rdname MultiFactor
#' @aliases names,anansi::MultiFactor-method
#' @export
#' @usage NULL
#'
S7::method(names, MultiFactor) <- function(x) {
    lapply(x@index, base::names)
}

#' @name MultiFactor
#' @rdname MultiFactor
#' @aliases dimnames
#' @export
#' @usage NULL
#'
S7::method(dimnames, MultiFactor) <- function(x) {
    dimnames(x@map)
}

#' @name MultiFactor
#' @importFrom methods show
#' @importMethodsFrom methods show
#' @aliases show,anansi::MultiFactor-method
#' @rdname MultiFactor
#' @usage NULL
#' @export
#'
S7::method(show, MultiFactor) <- function(object) {
    cat(
        "An ", class(object),
        ",\n    ", NCOL(object),
        " feature types across ",
        NROW(object),
        " edge lists.\n\n",
        sep = ""
    )
    Matrix::printSpMatrix(object@map)

    cat(
        "\nValues represent unique feature names in that edge list.\n\n",
        "Levels:\n\n",
        sep = ""
    )
    id_w <- max(nchar(colnames(object)))
    nm_w <- max(nchar(nlevels(object)))
    for (id in colnames(object)) {
        num_lvs <- length(levels(object)[[id]])
        cat(
            format(id, width = id_w),
            " : ",
            format(num_lvs, width = nm_w),
            " Levels: ",
            sep = ""
        )

        if (num_lvs > 4L) {
            cat(
                levels(object)[[id]][1],
                levels(object)[[id]][2],
                "...",
                levels(object)[[id]][num_lvs],
                "\n",
                sep = " "
            )
        } else {
            cat(levels(object)[[id]], "\n", sep = " ")
        }
    }
    invisible(NULL)
}

#' @name MultiFactor
#' @rdname MultiFactor
#' @aliases levels.anansi::MultiFactor
#' @export
#' @usage NULL
#'
S7::method(levels, MultiFactor) <- function(x) {
    x@levels
}

#' @name MultiFactor
#' @rdname MultiFactor
#' @aliases unfactor,anansi::MultiFactor-method
#' @export
#' @usage NULL
#'
S7::method(unfactor, MultiFactor) <- function(x) {
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
}

#' @name MultiFactor
#' @rdname MultiFactor
#' @description Analogous to `factors`. `droplevels(MultiFactor)` returns a
#'     `MultiFactor` with unused levels removed.
#' @param exclude `NULL` or `Named character list` of similar structure as
#'     `levels(MultiFactor)`. Which levels to drop from output.
#' @param select `NULL` or `Named character list` of similar structure as
#'     `levels(MultiFactor)`. Which levels to keep in output.
#' @details Only one of `select` and `exclude` should be provided, as they are
#'     each others complement.
#' @aliases droplevels.anansi::MultiFactor
#' @returns A MultiFactor
#' @export
#' @usage NULL
#'
S7::method(droplevels, MultiFactor) <- function(x, ..., exclude = NULL, select = NULL) {
    stopifnot(
        "Only one of 'exclude' and 'select' may be provided" = sum(
            is.null(exclude),
            is.null(select)
        ) >
            0L
    )
    stopifnot("'x' is not a MultiFactor." = is(x, "anansi::MultiFactor"))
    # Section 1. Trimming the indices by user input
    lvs <- levels(x)
    d <- x@map
    if (!is.null(exclude)) {
        stopifnot(
            "`'exclude' must be a named list of character vectors ." = is.list(
                exclude
            ) &&
                any(names(exclude) %in% names(lvs))
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
            "'select' arg must be a named list of character vectors." = is.list(
                select
            ) &&
                any(names(select) %in% names(lvs))
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
        # Loop over cols. First determine which rows are relevant per col/type
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
    return(x)    }

#' @name MultiFactor
#' @rdname MultiFactor
#' @export
#' @usage NULL
#'
S7::method(levels, MultiFactor) <- function(x) {
    x@levels
}
#' @name MultiFactor
#' @rdname MultiFactor
#' @param value a replacement character vector of suitable dimensions.
#' @export
#' @usage NULL
#'
S7::method(`levels<-`, MultiFactor) <- function(x, value) {
    x@levels <- value
    x
}

#' @name MultiFactor
#' @rdname MultiFactor
#' @param drop Whether to return a `list` (Default) or `MultiFactor`.
#' @export
#' @aliases [.anansi::MultiFactor
#' @usage NULL
#'
`[.anansi::MultiFactor` <- function(x, ..., drop = TRUE) {

    dot_args <- rlang::dots_list(
        ..., .preserve_empty = TRUE, .ignore_empty = "none"
    )
    dot_len <- length(dot_args)
    stopifnot("Too many arguments provided" = dot_len %in% seq(0L, 2L, 1L))
    missing_i <- rlang::is_missing(dot_args[[1L]])

    if(dot_len == 0L || (dot_len == 1L && missing_i)) { return(x) }

    d <- x@map
    l <- levels(x)
    x <- x@index

    if (dot_len == 1L) {
        ii <- rownames(d[dot_args[[1L]], , drop = FALSE])
        x <- x[ii]
    }
    if (dot_len == 2L) {
        missing_j <- rlang::is_missing(dot_args[[2L]])
        if (!missing_i) ii <- rownames(d[dot_args[[1L]], , drop = FALSE])
        if (!missing_j) jj <- colnames(d[, dot_args[[2L]], drop = FALSE])

        if (missing_i) {
            ii <- rowsWithCol(d, jj, FALSE)
            x <- lapply(x[ii], `[`, i = jj)
        } else if (missing_j) {
            x <- x[ii]
        } else {
            x <- lapply(x[ii], `[`, i = jj)
        }
    }

    if (drop) {
        return(x)
    }

    MultiFactor(x, levels = l)

}


#' @export
#' @name MultiFactor
#' @rdname MultiFactor
#' @aliases [<-.anansi::MultiFactor
#' @usage NULL
#'
`[<-.anansi::MultiFactor` <- function(
        x,
        ...,
        value
) {
    dot_args <- rlang::dots_list(
        ..., .preserve_empty = TRUE, .ignore_empty = "none"
    )
    dot_len <- length(dot_args)
    stopifnot("Too many arguments provided" = dot_len %in% seq(0L, 2L, 1L))

    if(dot_len == 0L) { return(x) }

    d <- x@map

    missing_i <- rlang::is_missing(dot_args[[1]])
    missing_j <- if (dot_len == 1L) TRUE else {
        rlang::is_missing(dot_args[[2L]])
    }

    if (!missing_i) ii <- rownames(d[dot_args[[1]], , drop = FALSE])
    if (!missing_j) jj <- colnames(d[, dot_args[[2]], drop = FALSE])

    if (missing_j) {
        x@index[ii] <- value
        return(x)
    }

    if (missing_i) {
        ii <- rowsWithCol(d, jj, names = TRUE)
    }

    for (i in ii) {
        for (j in jj) {
            x@index[[i]][, j] <- value[[i]][, j]
        }
    }
    return(x)
}

#' @export
#' @rdname MultiFactor
#' @name MultiFactor
#' @aliases `[[.anansi::MultiFactor`
#' @usage NULL
#'
`[[.anansi::MultiFactor` <- function(x, ...) {
    i <- rlang::dots_list(
        ..., .preserve_empty = TRUE, .ignore_empty = "none"
    )
    i_len <- length(i)
    stopifnot("Too many arguments provided" = i_len %in% seq(0L, 1L, 1L))

    # Empty returns self
    if(i_len == 0L) {return(x)}
    i <- i[[1L]]
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
}

#' @export
#' @rdname MultiFactor
#' @name MultiFactor
#' @aliases [[<-.anansi::MultiFactor
#' @usage NULL
#'
`[[<-.anansi::MultiFactor` <- function(x, ..., value) {
    i <- rlang::dots_list( ..., .preserve_empty = TRUE, .ignore_empty = "none")
    stopifnot("exactly one indexing value is required." = length(i) == 1L)
    i <- i[[1L]]
    d <- x@map
    # If i can't index d, stop. Appending not supported through `[[<-`.
    if (!all(i %in% colnames(d))) {
        if (anyNA(colnames(d)[i]) || length(colnames(d)[i]) != length(i)) {
            stop("No levels corresponding to `i` found in MultiFactor. ")
        }
    }
    ii <- rowsWithCol(d, i, FALSE)
    x@index[ii] <- value
    x
}


##############################################################################
##############################################################################
##############################################################################

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
    return(c(x, y[!names(y) %in% i]))
}


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
