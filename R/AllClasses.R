#' MultiFactor S7 container class
#' @name MultiFactor
#' @rdname MultiFactor-class
#' @description
#' `MultiFactor` is an S7 class to organize and manage multiple sets of factors,
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
#' @slot index Named `list` of named integer data frames of at least two columns
#'     each. The column names correspond to names in the `levels` slot. Similar
#'     to `factor`s, the integers in those columns correspond to the characters
#'     in that level. Accessed through regular list methods (e.g., `[`, `[[`).
#' @slot levels `Named list of character vectors`. Accessed through `levels(x)`
#' @slot map `(sparse) Matrix` specifying which elements contain which levels.
#'     Accesses through `x@dictionary`.
#' @returns a MultiFactor object.
#' @importFrom methods getClass
#' @inheritParams MultiFactor-methods
#' @seealso [MultiFactor-methods()]
#' @examples
#' # Generate some random linkage input
#' x <- data.frame(
#'     a = sample(letters[seq(3)], 10, replace = TRUE),
#'     A = sample(LETTERS[seq(3)], 10, replace = TRUE)
#' )
#' MultiFactor(x)
#'
#' @export
#'
MultiFactor <- S7::new_class(
    "MultiFactor",
    package = "anansi",
    properties = list(
        index = S7::class_list,
        levels = S7::class_list,
        map = methods::getClass("Matrix", where = "Matrix")
    ),
    constructor = function(x, levels = NULL, drop.unmatched = FALSE) {
        if (!is(x, "anansi::MultiFactor")) {
            if (validLinkDF(x)) {
                x <- list(x = x)
            }
            stopifnot(
                "Input not correctly formatted." = all(
                    vapply(
                        as.list(x, use.names = FALSE),
                        validLinkDF,
                        NA,
                        USE.NAMES = FALSE
                    )
                )
            )
        }
        if (is(x, "anansi::MultiFactor")) {
            if (is.null(levels)) {
                levels <- levels(x)
            }
            x <- x@index
        }
        x <- checkMergers(x)
        if (drop.unmatched) {
            x <- trimMultiFactor(x)
        }
        m <- mapMultiFactor(x)

        # Integer DF Input
        if (all(vapply(x, validIntLinkDF, NA, USE.NAMES = FALSE))) {
            stopifnot(
                "Input is integers, levels must be provided. " = !is.null(levels)
            )
            # Factor DF input
        } else if (all(vapply(x, validFactLinkDF, NA, USE.NAMES = FALSE))) {
            if (is.null(levels)) {
                levels <- factorInputMultiFactorLevels(x, m)
            }
            x <- listFactRefactor(x, m, levels)
            x <- lapply(x, factToIntDF)
            # Character DF input
        } else if (all(vapply(x, validCharLinkDF, NA, USE.NAMES = FALSE))) {
            if (is.null(levels)) {
                levels <- generateMultiFactorLevels(x, m)
            }
            x <- listCharToIntegers(x, m, levels)
        }
        # Get rid of row.names.
        x <- lapply(x, `row.names<-.data.frame`, value = NULL)
        S7::new_object(
            S7::S7_object(),
            index = x,
            levels = levels,
            map = m
        )
    }
)
S7::S4_register(MultiFactor)


#' @rdname MultiFactor-class
#' @name asMultiFactor
#' @param levels an optional named list of vectors of the unique values (as
#'     character strings) that x might have taken. The default is the unique set
#'     of values taken by lapply(x, as.character), sorted into increasing order
#'     of x.
#' @param drop.unmatched `Logical scalar` If `TRUE` (Default), for feature types
#'     that are seen at least twice, exclude features that only present in one
#'     of their respective link data frames.
#' @usage
#' ## Constructor for `MultiFactor` objects
#' MultiFactor(x, levels = NULL, drop.unmatched = FALSE)
#'
#' @export
#' @seealso \itemize{
#' \item [kegg_link()]: for an example of valid input.
#' }
#'
asMultiFactor <- function(x, levels = NULL, drop.unmatched = TRUE) {
    MultiFactor(x, levels, drop.unmatched)
}

#' @title AnansiWeb S7 container class
#' @rdname AnansiWeb-class
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
#' @slot tableY,tableX Two `matrix` objects of measurements, data. Rows are
#'     samples and columns are features. Access with `@tableY` and `@tableX`.
#' @slot dictionary `Matrix`, binary adjacency matrix. Optionally sparse.
#'     Typically generated using the `weaveWeb()` function. Access with
#'     `@dictionary`.
#' @slot metadata Optional `data.frame` of sample metadata. Access with
#'     `@metadata`.
#' @param tableY,tableX A table containing features of interest. Rows should be
#'     samples and columns should be features. Y and X refer to the position of
#'     the features in a formula: Y ~ X.
#' @param dictionary A binary adjacency matrix of class `Matrix`, or coercible
#'     to `Matrix`
#' @param metadata `list` of metadata. Optional.
#' @importClassesFrom Matrix Matrix
#' @importFrom Matrix Matrix drop0
#' @importFrom S4Vectors DataFrame
#' @export
#' @seealso \itemize{
#'  \item [AnansiWeb-methods]
#'  \item [randomAnansi]: For generation of random AnansiWeb objects.
#'  \item [kegg_link()]: For examples of input for link argument.
#' }
#' @returns an `AnansiWeb` object, with sparse binary biadjacency matrix
#' with features from `y` as rows and features from `x` as columns in
#' `dictionary` slot.
#' @usage
#' ## Constructor for `AnansiWeb` objects
#' AnansiWeb(tableX, tableY, dictionary, metadata = data.frame())
#' @examples
#'
#' # Use AnansiWeb() to construct an AnansiWeb object from components:
#' tX <- `dimnames<-`(replicate(5, (rnorm(36))),
#'     value = list(
#'         as.character(seq_len(36)),
#'         letters[1:5]
#'     )
#' )
#' tY <- `dimnames<-`(replicate(3, (rnorm(36))),
#'     value = list(
#'         as.character(seq_len(36)),
#'         LETTERS[1:3]
#'     )
#' )
#'
#' d <- matrix(TRUE,
#'     nrow = NCOL(tY), ncol = NCOL(tX),
#'
#'     # Note: Dictionary should have named dimensions
#'     dimnames = list(
#'         y_names = colnames(tY),
#'         x_names = colnames(tX)
#'     )
#' )
#' web <- AnansiWeb(tableX = tX, tableY = tY, dictionary = d)
#'
AnansiWeb <- S7::new_class(
    "AnansiWeb",
    package = "anansi",
    properties = list(
        tableY = S7::class_numeric | S7::class_logical,
        tableX = S7::class_numeric | S7::class_logical,
        dictionary = S7::class_any,
        metadata = S7::class_data.frame
    ),
    constructor = function(tableX, tableY, dictionary, metadata = data.frame()) {
        # coerce
        if (!is(dictionary, "Matrix")) {
            dictionary <- drop0(Matrix(dictionary, sparse = TRUE))
        }
        if (!is(tableX, "matrix")) tableX <- as.matrix(tableX)
        if (!is(tableY, "matrix")) tableY <- as.matrix(tableY)

        # check validity
        stopifnot(
            "'tableX' and 'tableY' need same number of rows (observations)" = NROW(
                tableX
            ) ==
                NROW(tableY)
        )
        stopifnot(
            "cols in 'tableY' need same amount as rows in dictionary" = NCOL(
                tableY
            ) ==
                NROW(dictionary)
        )
        stopifnot(
            "cols in 'tableX' need same amount as rows in dictionary" = NCOL(
                tableX
            ) ==
                NCOL(dictionary)
        )
        if (
            is.null(names(dimnames(dictionary))) ||
                any(names(dimnames(dictionary)) %in% "")
        ) {
            warning("Dimnames of 'dictionary' were missing; Assigned 'y' and 'x'.")
            names(dimnames(dictionary)) <- c("y", "x")
        }
        metadata <- .check_metadata_labels( metadata, tableY, tableX )
        print(metadata)


        # return AnansiWeb
        S7::new_object(
            S7::S7_object(),
            tableY = `rownames<-`(tableY, row.names(metadata)),
            tableX = `rownames<-`(tableX, row.names(metadata)),
            dictionary = dictionary,
            metadata = as.data.frame(metadata)
        )
    }
)
S7::S4_register(AnansiWeb)

#' AnansiWeb constructor helper function
#' @noRd
.check_metadata_labels <- function(metadata, tableY, tableX) {
    ry <- if(is.null(row.names(tableY))) {
        as.character(seq_len(NROW(tableY))) } else row.names(tableY)
    rx <- if(is.null(row.names(tableX))) {
        as.character(seq_len(NROW(tableX))) } else row.names(tableX)

    id <- if(identical(ry, rx)) ry else as.character(seq_along(ry))

    rn <- id
    prepend_these <- !grepl("^anansi_ID_", rn)
    rn[prepend_these] <- paste0("anansi_ID_", rn[prepend_these])

    if(rlang::is_empty(metadata)) {
        warning(
            "Argument `metadata` not provided; Please validate sample ID order."
        )
        metadata <- data.frame(row.names = rn)
        return(metadata)
    }

    if(NROW(metadata) != length(rn)) {
        stop("Metadata must have exactly one row per sample in tableY,tableX")
        }

    row.names(metadata) <- rn

    return(metadata)
}

#' @title AnansiTale S7 container class. Not intended for general use.
#' @rdname AnansiTale
#' @name AnansiTale
#' @aliases AnansiTale-class
#' @slot subject A character that describes the data that was queried.
#' @slot type A character that describes type of parameter contained in the
#'     `estimates` slot. For example r.values for correlations or r.squared
#'     for models.
#' @slot df a vector of length 2, containing df1 and df2 corresponding to the
#'     F-ratio considered.
#' @slot estimates A matrix containing the estimates for the parameters named in
#'     the `type` slot.
#' @slot f.values A matrix containing the f-values, for least-squares.
#' @slot t.values A matrix containing the t-values, for correlations.
#' @slot p.values A matrix containing the p.values for the parameters named in
#'     the `type` slot.
#' @description `AnansiTale` is the main container that will hold your
#'     stats output data coming out of the `anansi` pipeline.
#' @param subject A character that describes the data that was queried.
#' @param type A character that describes type of parameter contained in the
#'     `estimates` slot. For example r.values for correlations or r.squared
#'     for models.
#' @param df a vector of length 2, containing df1 and df2 corresponding to the
#'     F-ratio considered.
#' @param estimates A matrix containing the estimates for the parameters named in
#'     the `type` slot.
#' @param f.values A matrix containing the f-values, for least-squares.
#' @param t.values A matrix containing the t-values, for correlations.
#' @param p.values A matrix containing the p.values for the parameters named in
#'     the `type` slot.
#' @returns an `AnansiTale`
#' @examples
#' AnansiTale
#'
#' @export
#'
AnansiTale <- S7::new_class(
    "AnansiTale",
    package = "anansi",
    properties = list(
        subject = S7::class_character,
        type = S7::class_character,
        df = S7::class_numeric,
        estimates = S7::class_numeric,
        f.values = S7::class_numeric,
        t.values = S7::class_numeric,
        p.values = S7::class_numeric
    )
)
