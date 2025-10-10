#' Methods for MultiFactor S7 container class
#' @name MultiFactor-methods
#' @rdname MultiFactor-methods
#' @aliases levels.anansi::MultiFactor
#' @aliases unfactor unfactor,anansi::MultiFactor-method
#' @description \describe{
#'     \item{`droplevels()` }{`droplevels(MultiFactor)` returns a `MultiFactor`
#'     with unused levels removed, Analogous to the `factor` method.}
#'     }
#' @examples
#' # Setup
#' x <- MultiFactor(kegg_link())
#' x
#'
#' # Basic properties
#' dim(x)
#' dimnames(x)
#'
#' # Factor-like properties
#' head(levels(x)$ko)
#' droplevels(x)
#' head(unfactor(x)$ec2ko)
#'
#' # Extract common output formats
#' getEdgeList(x)
#'
#' @param x `MultiFactor` on which the method should be applied, or, in
#'     case of the constructor `MultiFactor()`, a named `list` of data.frames
#'     with two named columns each, where elements that share a row indicates
#'     thet are adjacent.
#' @param ... `i,j` indices specifying elements to extract or replace. Indices
#'     are numeric or character vectors or empty (missing) or NULL. Numeric
#'     values are coerced to integer or whole numbers as by as.integer or for
#'     large values by trunc (and hence truncated towards zero). Character
#'     vectors will be matched to the names of the object.
#' @param value Replacement value, typically of same type as that which is to be
#'     replaced.
#' @param exclude `NULL` or `Named character list` of similar structure as
#'     `levels(MultiFactor)`. Which levels to drop from output.
#' @param select `NULL` or `Named character list` of similar structure as
#'     `levels(MultiFactor)`. Which levels to keep in output.
#' @param use.names,ignore.mcols For compatibility, not used.
#' @details Only one of `select` and `exclude` should be provided, as they are
#'     each others complement.
#' @examples
#' droplevels(x, exclude = list(ko = "K00001"))
#' droplevels(x, select = list(ko = "K00001"))
#' @returns A MultiFactor
NULL

#' @export
#'
S7::method(getEdgeList, MultiFactor) <- function(x) {
    as.data.frame(do.call(rbind, base::names(x)))
}

#' @export
#'
S7::method(dim, MultiFactor) <- function(x) dim(x@map)


#' @export
#'
S7::method(names, MultiFactor) <- function(x) {
    lapply(x@index, base::names)
}

#' @export
#'
S7::method(dimnames, MultiFactor) <- function(x) {
    dimnames(x@map)
}

#' @importMethodsFrom methods show
#' @export
#'
S7::method(show, MultiFactor) <- function(object) {
    cat(
        "An ", paste(class(object), collapse = " "),
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

#' @export
#'
S7::method(print, MultiFactor) <- function(x, ...) show(x)


#' @importMethodsFrom S4Vectors unfactor
#' @importFrom S4Vectors unfactor
#' @export
#'
S7::method(unfactor, MultiFactor) <-
    function(x, use.names = TRUE, ignore.mcols = TRUE) {
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

#' @export
#' @importMethodsFrom S4Vectors droplevels
#'
S7::method(droplevels, MultiFactor) <- function(
        x, ..., exclude = NULL, select = NULL
) {
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
    return(x)
}


#' @export
S7::method(levels, MultiFactor) <- function(x) `levels.anansi::MultiFactor`(x)

#' @export
`levels.anansi::MultiFactor` <- function(x) {
    x@levels
}

local({
    S7::method(`[`, MultiFactor) <- function(x, i, j, ..., drop = TRUE) {
        if (!all(names(sys.call()) %in% c("", "drop"))) {
            warning("named arguments other than 'drop' are discouraged")
        }
        raw_call <- rlang::call_match(
            dots_expand = FALSE, defaults = TRUE
        )


        do.call(
            ".sub_MultiFactor", rlang::call_args(raw_call),
            quote = FALSE, envir = rlang::caller_env()
        )
    }


    S7::method(`[[`, MultiFactor) <- function(x, i, ...) {
        raw_call <- rlang::call_match(
            dots_expand = FALSE, defaults = TRUE, fn = .subsub_MultiFactor
        )

        do.call(".subsub_MultiFactor", rlang::call_args(raw_call),
                envir = rlang::caller_env()
        )
    }
})

#' @export
`[<-.anansi::MultiFactor` <- function(x, i, j, ..., value) {
    if (!all(names(sys.call()) %in% c("", "value"))) {
        warning("named arguments are discouraged")
    }

    raw_call <- rlang::call_match(
        dots_expand = FALSE, defaults = TRUE, fn = .sub_rep_MultiFactor
    )
    do.call(".sub_rep_MultiFactor", rlang::call_args(raw_call),
            envir = rlang::caller_env()
    )
}

#' @export
`[[<-.anansi::MultiFactor` <- function(x, i, ..., value) {
    if (!all(names(sys.call()) %in% c("", "value"))) {
        warning("named arguments are discouraged")
    }
    raw_call <- rlang::call_match(
        dots_expand = FALSE, defaults = TRUE, fn = .sub_sub_rep_MultiFactor
    )
    do.call(".sub_sub_rep_MultiFactor", rlang::call_args(raw_call),
            envir = rlang::caller_env()
    )
}



.sub_MultiFactor <- function(x, i, j, drop = TRUE) {
    missing_i <- rlang::is_missing(i)
    missing_j <- rlang::is_missing(j)
    missing_x <- rlang::is_missing(x)
    dot_len <- sum(!missing_i, !missing_j)

    stopifnot("Too many arguments provided" = dot_len %in% seq(0L, 2L, 1L))
    # missing_i <- rlang::is_missing(dot_args[[1L]])

    if (dot_len == 0L) {
        return(x)
    }

    d <- x@map
    l <- levels(x)
    x <- x@index

    if (!missing_i) {
        ii <- rownames(d[i, , drop = FALSE])
        x <- x[ii]
    }
    if (dot_len == 2L) {
        if (!missing_i) ii <- rownames(d[i, , drop = FALSE])
        if (!missing_j) jj <- colnames(d[, j, drop = FALSE])

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


#' @importFrom rlang is_missing
.subsub_MultiFactor <- function(x, i) {
    # Empty returns self
    if (rlang::is_missing(i)) {
        return(x)
    }
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
    .sub_MultiFactor(x, ii)
}


.sub_rep_MultiFactor <- function(x, i, j, value) {
    missing_i <- rlang::is_missing(i)
    missing_j <- rlang::is_missing(j)
    missing_x <- rlang::is_missing(x)
    dot_len <- sum(!missing_i, !missing_j)

    stopifnot("Too many arguments provided" = dot_len %in% seq(0L, 2L, 1L))

    if (dot_len == 0L) {
        return(x)
    }

    d <- x@map

    if (!missing_i) ii <- rownames(d[i, , drop = FALSE])
    if (!missing_j) jj <- colnames(d[, j, drop = FALSE])

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

.sub_sub_rep_MultiFactor <- function(x, i, value) {
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

#' For S4 compatibility
#' @noRd
unfactor <- function(x) S4Vectors::unfactor(x)

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
