#' Weave an AnansiWeb
#' @rdname weaveWeb
#' @order 0
#' @param formula `formula` of the form y ~ x, denoting desired output
#'     format; assigns y to rows and columns to x. Equivalent to using `x`
#'     and `y` arguments.
#' @param x,y `Character scalar`, names of feature types that should be
#'     linked. Should be found in the column names of `link`.
#' @param link One of the following:
#' \itemize{
#'  \item `Character scalar` with value `"none"`.
#'  \item `data.frame` with two columns
#'  \item `list` with two such `data.frame`s.
#' }
#' @param tableY,tableX A table containing features of interest. Rows should be
#'     samples and columns should be features. Y and X refer to the position of
#'     the features in a formula: Y ~ X.
#' @param ... further arguments.
#' @details
#' If the `link` argument is `"none"`, all features will be considered
#' linked. If one or more `data.frame`s, colnames should be as specified in
#' `x` and `y`.
#' @seealso \itemize{
#'  \item [AnansiWeb-methods]: For utility functions to get and set.
#'  \item [AnansiWeb()]: For more general constructor.
#'  \item [kegg_link()]: For examples of input for link argument.
#'  \item [getWeb()]: For
#'  [MultiAssayExperiment::MultiAssayExperiment()] methods.
#' }
#'
#' @returns an `AnansiWeb` object, with sparse binary biadjacency matrix
#' with features from `y` as rows and features from `x` as columns in
#' `dictionary` slot.
#' @description
#' Generate a biadjacency matrix, linking the features between two tables.
#' Return an `AnansiWeb` object which contains all three.
#'
#' `weaveWeb()` is for general use and has flexible default settings.
#'
#' `weaveKEGG()` is a wrapper that sets `link` to `kegg_link()`.
#' All variants are special cases of `weaveWeb()`.
#' @examples
#' # Basic usage
#' weaveWeb(cpd ~ ko, link = kegg_link())
#' weaveWeb(x = "ko", y = "ec", link = ec2ko)
#' weaveWeb(ec ~ cpd, link = ec2cpd)
#'
#' # A wrapper is available for kegg ko, ec and cpd data
#' generic      <- weaveWeb(cpd ~ ko, link = kegg_link())
#' kegg_wrapper <- weaveKEGG( cpd ~ ko )
#'
#' identical(generic, kegg_wrapper)
#'
#' # The following are equivalent to transposition:
#' a <- weaveWeb(ko ~ cpd, link = kegg_link())$dictionary
#' b <- weaveWeb(cpd ~ ko, link = kegg_link())$dictionary
#'
#' identical(a, Matrix::t(b))
#'
#' @usage NULL
#' @export
#'
weaveWeb <- function(x, ...) UseMethod("weaveWeb")

#' @rdname weaveWeb
#' @importFrom Matrix Matrix
#' @order 2
#' @export
#'
weaveWeb.default <- function(x, y, link = NULL, tableX = NULL, tableY = NULL,
                             metadata = NULL, ...){
    terms <- c(y, x)
    stopifnot("both 'x' and 'y' terms must be provided as character" =
                  is(terms, "character") && length(terms) == 2L)
    if(identical(link, "none")) return(web_missing_link(tableX, tableY, terms))

    # Ensure link is a MultiFactor
    link  <- MultiFactor(link)

    stopifnot("both 'x' and 'y' terms must be found as colnames in 'link'" =
                  all(terms %in% colnames(link)) )
    # Determine required ids in order
    all_terms <- termSeq(x, y, link)

    # Trim link levels and tables based on feature overlap
    if(!is.null(tableX)) {
        link   <- trimByInput(link, tableX, x)
        tableX <- tableX[, sort(intersect(colnames(tableX), levels(link)[[x]]))]
    }
    if(!is.null(tableY)) {
        link   <- trimByInput(link, tableY, y)
        tableY <- tableY[, sort(intersect(colnames(tableY), levels(link)[[y]]))]
    }
    # Construct dictionary
    d <- dictionaryMatrix(link, all_terms)
    dimnames(d) <- list(y = colnames(tableY), x = colnames(tableX))
    names(dimnames(d)) <- c(y, x)

    # Dummy tables if missing
    if(is.null(tableX) && is.null(tableY)) {
        dimnames(d) <- levels(link)[terms]
        tableY = matrix(ncol = NROW(d), dimnames = list(NULL, rownames(d)))
        tableX = matrix(ncol = NCOL(d), dimnames = list(NULL, colnames(d)))
    }
    #
    AnansiWeb(
        tableY     = as.matrix(tableY)[,rownames(d), drop = FALSE],
        tableX     = as.matrix(tableX)[,colnames(d), drop = FALSE],
        dictionary = d,
        metadata   = metadata)
}

#' @rdname weaveWeb
#' @export
#' @order 1
#'
weaveWeb.formula <- function(
        formula, link = NULL, tableX = NULL, tableY = NULL, ...
) {
    if (missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")

    terms <- all.vars(formula)
    if(is.null(link) || identical(link, "none")) return(
        weaveWeb.default(x = terms[2], y = terms[1], link, tableX, tableY)
    )

    link <- MultiFactor(link)

    if (sum(terms %in% colnames(link)) != 2L)
        stop("Variables from 'formula' not found in 'link'.")

    weaveWeb.default(x = terms[2], y = terms[1], link, tableX, tableY)
}

#' @rdname weaveWeb
#' @export
#'
weaveKEGG <- function(x, ...) weaveWeb(x, link = kegg_link(), ...)



###############################################################################
###############################################################################

#' Find a path through different feature types.
#' @inheritParams deliver_web_single
#' @importFrom igraph shortest_paths
#' @returns a Character vector of the ids to walk in order.
#' @noRd
#'
termSeq <- function(x, y, link){
    g <- getGraph(link)
    sp <- shortest_paths(g, from = y, to = x, output = "vpath")
    names(unlist(sp, FALSE, FALSE)[[1]])
}

#' Find the order in which link data frames should be listed
#' @param term_list list of `Character vectors`, each with length of two.
#' @param link_names names(link)
#' @returns a numeric vector with order in which row data frames should be
#'     traversed.
#' @noRd
#'
stepSeq <- function(term_list, link_names){
    vapply(term_list,
           function(st) which(
               vapply(link_names,
                      function(x) all(st %in% x),
                      FUN.VALUE = TRUE, USE.NAMES = FALSE)),
           FUN.VALUE = 0L, USE.NAMES = FALSE)
}

#' Generate dictionary Matrix from link input
#' @inheritParams weaveWeb
#' @importMethodsFrom Matrix %&%
#' @noRd
#'
dictionaryMatrix <- function(link, all_terms){
    term_list <- lapply(seq_len(length(all_terms)-1L),
                        FUN = function(x) all_terms[c(x, x + 1L)])
    lv_len    <- vapply(levels(link), length, 0L, USE.NAMES = TRUE)
    lv_list   <- lapply(term_list, function(x) lv_len[x])
    steps     <- stepSeq(term_list, names(link))

    # Handle simple case of one link df first, return sparse matrix.
    if(length(steps) == 1L)
        return(
            mapFromLink(all_terms, df = link[[steps]], dims = lv_len[all_terms])
            )

    # Otherwise, make a list of matrices to Reduce to final dictionary
    mat_list <- mapply(mapFromLink,
                       terms = term_list,
                       df = link[steps],
                       dims = lv_list)
    Reduce(`%&%`,  mat_list)

}

#' @param terms y, x id of cols
#' @param df element of a `MultiFactor` object
#' @param dims length-2 integer vector of matrix dimensions.
#' @importFrom Matrix sparseMatrix
#' @returns a sparse biadjacency Matrix
#' @noRd
#'
mapFromLink <- function(terms, df, dims)
    sparseMatrix(i = df[[terms[1]]], j = df[[terms[2]]], dims = dims)

#' @description Called by weaveWeb to subset link to inly include the features
#'     found in the input table.
#' @returns a MultiFactor subsetted by relevant features
#' @param link a `MultiFactor` .
#' @param id `Character scalar`, naming the x term to be trimmed
#' @param tableID A table containing features of interest, `tableX` or `tableY`.
#' @noRd
#'
trimByInput <- function(link, tableID, id) {
    x.names <- names(link)

    x.ind <- vapply(x.names, `%in%`, x = id, NA, USE.NAMES = FALSE)
    sel.obj <- link[x.ind]
    term_list <- lapply(sel.obj, function(df) df[, id])
    term_list[["table_IDs"]] <- match( colnames(tableID), levels(link)[[id]] )
    keep    <- Reduce(intersect, term_list)

    # Filter feature ids in each df to only include universally shared ones.
    link[x.ind] <- lapply(sel.obj, function(df) {
        df <-  df[df[[id]] %in% keep,]
        df[,id] <- as.integer(factor(df[,id]))
        return(df)})
    levels(link)[[id]] <- levels(link)[[id]][keep]

    link
}


#' Produce a biadjacency matrix given tables and a dictionary
#' @description calculates a biadjacency matrix for the cases where
#' `link` is a single - or a `list` of two - `data.frame`(s).
#' @param link a `data.frame` or `list` of two compatible ones.
#' @param terms a length 2 character vector, naming the x and y terms in order
#' @param tableY A table containing features of interest. Rows should be samples
#' and columns should be features. The Y and X refer to the position of the
#' features in a formula: Y ~ X.
#' @param tableX A table containing features of interest. Rows should be samples
#' and columns should be features. The Y and X refer to the position of the
#' features in a formula: Y ~ X.
#' @returns a sparse boolean biadjacency matrix, for use in main anansi workflow
#' @noRd
#'
deliver_web_single <- function(link, terms, tableX, tableY){

    lv <- levels(link)[terms]
    link <- link[[1]]
    for(id in seq_along(terms))  {
        link[[terms[id]]] <-  lv[[terms[id]]][link[[terms[id]]]]
    }

    if(!is.null(tableX))  link <- link[link[,1L] %in% colnames(tableX),]
    if(!is.null(tableY))  link <- link[link[,2L] %in% colnames(tableY),]

    d <- df_to_sparse_biadjacency_matrix(link)
    names(dimnames(d)) <- rev(terms)

    return( d )

}

#' Produce a biadjacency matrix given tables and a dictionary
#' @description calculates a biadjacency matrix for the cases where
#' `link` is a single - or a `list` of two - `data.frame`(s).
#' @param link a `data.frame` or `list` of two compatible ones.
#' @param terms a length 2 character vector, naming the x and y terms in order
#' @param tableY A table containing features of interest. Rows should be samples
#' and columns should be features. The Y and X refer to the position of the
#' features in a formula: Y ~ X.
#' @param tableX A table containing features of interest. Rows should be samples
#' and columns should be features. The Y and X refer to the position of the
#' features in a formula: Y ~ X.
#' @returns a sparse boolean biadjacency matrix, for use in main anansi workflow
#' @noRd
#'
deliver_web_list <- function(link, terms, tableX, tableY){


    # order row_col:
    # 1_2, 2_3, 3_4, 4_5, 5_6 ... N-1_N


    #
    #
    # A_B <- matrix(T, ncol = 20, nrow = 40, dimnames = list("A" = paste("a", 1:40, sep = "_"),
    #                                                        "B" = paste("b", 1:20, sep = "_")))
    # B_C <- matrix(T, ncol = 30, nrow = 20, dimnames = list("B" = paste("b", 1:20, sep = "_"),
    #                                                        "C" = paste("c", 1:30, sep = "_")))
    # C_D <- matrix(T, ncol = 50, nrow = 30, dimnames = list("C" = paste("c", 1:30, sep = "_"),
    #                                                        "D" = paste("d", 1:50, sep = "_")))
    #
    # mat_list <- list(A_B, B_C, C_D)
    #
    #
    # # This works!
    # Reduce(`%*%`, accumulate = TRUE,  mat_list)
    #


    cn.1 <- colnames(link[[1L]]); cn.2 <- colnames(link[[2L]])
    i <- intersect(cn.1,cn.2)
    stopifnot("data.frames in 'link' must share a colname" = length(i) == 1L)

    if(all(cn.1 %in% c(i, terms[1]))) {

        x.df <- link[[1L]][, c(i, terms[1L])]
        y.df <- link[[2L]][, c(i, terms[2L])]

    } else   if(all(cn.1 %in% c(i, terms[2]))) {

        x.df <- link[[2L]][, c(i, terms[1L])]
        y.df <- link[[1L]][, c(i, terms[2L])]

    }

    d <- web_from_2_dfs(x.df, y.df, colnames(tableX), colnames(tableY))
    names(dimnames(d)) <- rev(terms)
    return( d )
}

#' @importFrom Matrix tcrossprod
#' @noRd
#'
web_from_2_dfs <- function(x.df, y.df, x.ids, y.ids){
    if(!is.null(x.ids)) { x.df <- x.df[x.df[, 2L] %in% x.ids, ] }
    if(!is.null(y.ids)) { y.df <- y.df[y.df[, 2L] %in% y.ids, ] }

    i <- sort(intersect(x.df[,1L],y.df[,1L]))
    x.df <- x.df[x.df[,1L] %in% i,]
    y.df <- y.df[y.df[,1L] %in% i,]

    x.mat <- df_to_sparse_biadjacency_matrix(x.df)
    y.mat <- df_to_sparse_biadjacency_matrix(y.df)

    tcrossprod(y.mat, x.mat, boolArith = TRUE)
}

#' @importFrom igraph vertex_attr<- vertex_attr graph_from_data_frame as_biadjacency_matrix
#' @importFrom Matrix drop0
#' @noRd
#'
df_to_sparse_biadjacency_matrix <- function(x){
    x.g <- graph_from_data_frame(x, directed = FALSE)
    vertex_attr(x.g, name = "type") <- vertex_attr(x.g, "name") %in% x[,1L]
    m <- drop0(as_biadjacency_matrix(x.g, sparse = TRUE))
    m <- m[order(rownames(m)),order(colnames(m))]
    return(m)
}

#' Make a full web; for all vs all association testing
#' @description
#' Make a fully TRUE biadjacency matrix with dimensions of the two input tables.
#' @param tableX,tableY `matrix` of features of table `X`.
#' @param terms `character vector` names of x & y terms
#' @returns
#' An `AnansiWeb` object with both tables and a fully `TRUE`
#' (non-sparse) matrix from the `Matrix` package.
#' @importFrom Matrix Matrix
#' @noRd
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
        tableY     = as.matrix(tableY)[,rownames(d)],
        tableX     = as.matrix(tableX)[,colnames(d)],
        dictionary = d,
        metadata = metadata)

}


