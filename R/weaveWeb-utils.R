# Find a path through different feature types.
# returns a Character vector of the ids to walk in order.
#' @importFrom igraph shortest_paths
#'
termSeq <- function(x, y, link) {
    stopifnot(
        "both 'x' and 'y' terms must be found as colnames in 'link'" = all(
            c(x, y) %in% colnames(link)
        )
    )
    g <- getGraph(link)
    sp <- igraph::shortest_paths(g, from = y, to = x, output = "vpath")
    names(unlist(sp, FALSE, FALSE)[[1]])
}

#' Find the order in which link data frames should be listed
#' @param term_list list of `Character vectors`, each with length of two.
#' @param d link@map.
#' @returns a numeric vector with order in which row data frames should be
#'     traversed.
#' @noRd
#'
stepSeq <- function(term_list, d) {
    vapply(
        term_list,
        rowsWithCol,
        d = d,
        name = FALSE,
        FUN.VALUE = 0L,
        USE.NAMES = FALSE
    )
}

subsetByPath <- function(link, all_terms) {
    term_list <- lapply(
        seq_len(length(all_terms) - 1L),
        FUN = function(x) all_terms[c(x, x + 1L)]
    )
    steps <- stepSeq(term_list, link@map)
    link@index <- link[steps]
    link@levels <- link@levels[all_terms]
    link@map <- mapMultiFactor(link[steps], mode = "counts")

    return(link)
}

#' Generate dictionary Matrix from link input
#' @inheritParams weaveWeb
#' @param all_terms `Character vector` of all path terms in sequence.
#'     `termSeq(x, y, link)`
#' @importMethodsFrom Matrix %&%
#' @noRd
#'
dictionaryMatrix <- function(link, all_terms) {
    term_list <- lapply(
        seq_len(length(all_terms) - 1L),
        FUN = function(x) all_terms[c(x, x + 1L)]
    )
    steps <- stepSeq(term_list, link@map)
    lv_len <- vapply(X = levels(link), FUN = length, 0L, USE.NAMES = TRUE)

    # Handle simple case of one link df first, return sparse matrix.
    if (length(steps) == 1L) {
        return(mapFromLink(
            all_terms,
            df = link@index[[steps]],
            dims = lv_len[all_terms]
        ))
    }
    lv_list <- lapply(term_list, function(x) lv_len[x])
    # Otherwise, make a list of matrices to Reduce to final dictionary
    mat_list <- mapply(
        FUN = mapFromLink,
        terms = term_list,
        df = link@index[steps],
        dims = lv_list
    )
    Reduce(Matrix::`%&%`, mat_list)
}

#' @param terms id of cols. `c(y, x)`.
#' @param df element of a `MultiFactor` object
#' @param dims length-2 integer vector of matrix dimensions.
#' @importFrom Matrix sparseMatrix
#' @returns a sparse biadjacency Matrix
#' @noRd
#'
mapFromLink <- function(terms, df, dims) {
    Matrix::sparseMatrix(i = df[[terms[1]]], j = df[[terms[2]]], dims = dims)
}

#' @description Called by weaveWeb to subset link to inly include the features
#'     found in the input table.
#' @returns a MultiFactor subsetted by relevant features
#' @param link a `MultiFactor` .
#' @param id `Character scalar`, naming the x term to be trimmed
#' @param tableID A table containing features of interest, `tableX` or `tableY`.
#' @noRd
#'
trimByInput <- function(link, tableID, id) {
    lv <- levels(link)[[id]]
    d <- link@map
    r <- rowsWithCol(d, id)
    stopifnot(
        "Feature names appeared in several index elements. " = length(r) == 1L
    )
    # Subset index by table columns
    xr <- link@index[[r]]
    x.names <- match(colnames(tableID), lv)
    xr <- xr[xr[, id] %in% x.names, ]
    xr.id <- xr[, id]

    # Subset levels
    link@levels[[id]] <- lv[sort(unique(xr.id))]
    # Reorder and replace indices
    xr[, id] <- match(xr.id, sort(unique(xr.id)))
    link@index[[r]] <- xr

    link@map <- mapMultiFactor(link@index, mode = "counts")

    return(link)
}


# Make a full web; for all vs all association testing
# Make a fully TRUE biadjacency matrix with dimensions of the two input tables.
# tableX, tableY: matrix of features of tableX.
# terms: character vector names of x & y terms
# Returns an `AnansiWeb` object with both tables and a fully `TRUE`
# (non-sparse) matrix from the `Matrix` package.
#' @importFrom Matrix Matrix
#'
web_missing_link <- function(tableX, tableY, terms, metadata = NULL) {
    d <- Matrix(
        data = TRUE,
        nrow = NCOL(tableY),
        ncol = NCOL(tableX),
        dimnames = list(sort(colnames(tableY)), sort(colnames(tableX)))
    )
    names(dimnames(d)) <- rev(terms)

    AnansiWeb(
        tableY = as.matrix(tableY)[, rownames(d)],
        tableX = as.matrix(tableX)[, colnames(d)],
        dictionary = d,
        metadata = metadata
    )
}

# TRUE if i can select in x
# i: Character or numeric scalar. Index to check.
# x: object to check i in
# returns TRUE if i selects in x, FALSE otherwise
#
.check_valid_selection <- function(i, x) {
    # Need to be length 1.
    if (length(i) != 1L) FALSE

    # If numeric, needs to be within element length of x
    if (is.numeric(i)) i <= length(x)
    # If character, x needs to be named and i needs to be within those names
    if (is.character(i)) {
        if (is.null(names(x))) FALSE

        !is.na(match(i, names(x)))
    }
    # If that didn't work, invalid selection. return FALSE.
    FALSE
}

.check_fixed_args <- function(kwargs) {
    # Check fixed arguments
    fixed_args <- c("web", "metadata")
    remove <- names(kwargs) %in% fixed_args
    # If fixed arguments in kwargs, remove them
    if (any(remove)) {
        removed <- paste0(
            names(kwargs[remove]),
            sep = "'",
            collapse = ", '"
        )
        warning(
            "The arguments '",
            removed,
            " should not be used, ",
            "as they are extracted from 'x'.",
            call. = FALSE
        )
    }
    kwargs <- kwargs[!remove]
    return(kwargs)
}

#' @noRd
.test_coherent <- function(tab, exp, tab.type, ass.type) {
    e_out <- unique(c(tab, exp))
    if (length(e_out) == 0L) {
        e_out <- 1L
    }
    stopifnot(
        "args 'tableY,X' cannot contradict args 'experiment1,2'." =
            length( e_out ) == 1L
    )

    t_out <- unique(c(tab.type, ass.type))
    if (length(t_out) == 0L) {
        t_out <- 1L
    }
    stopifnot(
        "args 'typeY,X.' cannot contradict args 'assay.type1,2'." =
            length( t_out ) == 1L
    )

    return(list(e_out, t_out))
}

#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom SingleCellExperiment altExp altExpNames
.get_table_from_tse <- function(tse, y_id, x_id){
    all_assays <- c(assayNames(tse), altExpNames(tse))
    if(is.numeric(y_id)) { y_id <- all_assays[y_id] }
    if(is.numeric(x_id)) { x_id <- all_assays[x_id] }


    tab_list <- lapply(c(y_id, x_id),
           FUN = function(id){

               if( !id %in% all_assays ) {
                   stop("no assay with name ",id, " found.", call. = FALSE)
               }
               if( id %in% assayNames(tse) ){
                   return( tse )
               } else {
                   return( altExp(tse, id) )
               }
           }
    )
    names(tab_list) <- c(y_id, x_id)

    return(tab_list)
}

#
.test_mae_has_exp <- function (x, experiment) {
    if (is.numeric(experiment) && experiment > length(experiments(x))) {
        stop("'", deparse(substitute(experiment)), "' is greater than the ",
             "number of experiments in MAE object.", call. = FALSE)
    }
    if (!(is.character(experiment) && experiment %in% names(experiments(x)) ||
          is.numeric(experiment) && experiment <= length(experiments(x)))) {
        stop("'", deparse(substitute(experiment)), "' ", "must be numeric or character value specifying ",
             "experiment in experiment(x).", call. = FALSE)
    }
    obj <- x[[experiment]]
    if (!(is(obj, "SummarizedExperiment"))) {
        stop("The class of experiment specified by ", deparse(substitute(experiment)),
             " must be 'SummarizedExperiment'.", call. = FALSE)
    }
    return(NULL)
}
