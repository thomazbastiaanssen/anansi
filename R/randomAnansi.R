#' Generate a random AnansiWeb or AnansiLinkMap
#' @name randomAnansi
#' @description
#' Randomly generate a valid `AnansiWeb` or `AnansiLinkMap` object.
#' @param n_samples `Numeric scalar` Number of samples to be generated.
#' @param n_features_y,n_features_x `Numeric scalar` Number of features to be
#'     generated.
#' @param tableY,tableX A table containing features of interest. Rows should be
#'     samples and columns should be features. Y and X refer to the position of
#'     the features in a formula: Y ~ X.
#' @param dictionary A binary adjacency matrix of class `Matrix`, or
#' coercible to `Matrix`.
#' @returns a randomly generated object of the specified class.
#' @examples
#' # Make a random AnansiWeb object
#' randomWeb()
#' randomLinkMap()
#' @seealso [AnansiWeb()], [AnansiLinkMap()]
#'
NULL

#' @rdname randomAnansi
#' @name randomWeb
#' @export
#'
randomWeb <- function(n_samples = 10, n_features_x = 8, n_features_y = 12,
                      sparseness = 0.5, tableY = NULL, tableX = NULL,
                      dictionary = NULL){
    stopifnot("'sparseness' must be a proportion [0-1]. " =
                  sparseness <= 1 && sparseness > 0)
    stopifnot("At least one of 'tableY,tableX', 'dictionary' should be NULL." =
                  any(c(is.null(tableY), is.null(tableX), is.null(dictionary))))
    stopifnot("'tableY,tableX' should either both be provided or both NULL. " =
                  is.null(tableY) == is.null(tableX) )

    density <- 1 - sparseness
    # All missing: return full random Web
    if(all(c(is.null(tableY), is.null(tableX), is.null(dictionary)))) return(
        randomWebFull(n_samples, n_features_x, n_features_y, density) )
    # Dictionary missing: make random fitting dictionary, return filled Web
    if(is.null(dictionary)) return(
        randomWebDic(tableY, tableX, density) )
    # Tables missing: make random fitting tables, return filled Web
    return(randomWebTab(n_samples, dictionary))
}

#' @rdname randomAnansi
#' @param n_types `Numeric scalar`, number of types of features to generate
#' @param n_features `Numeric scalar`, number of features per type
#' @param sparseness `Numeric scalar`, proportion: How rare are connections
#' @export
#'
randomLinkMap <- function(n_types = 6, n_features = 100,
                          sparseness = 0.5){
    stopifnot("'sparseness' must be a proportion [0-1]. " =
                  sparseness <= 1 && sparseness > 0)
    n_types <- max(min(n_types, 26), 2)
    ids <- letters[seq_len(n_types)]
    out_names <- paste0(ids[-n_types], "2",ids[-1L])
    id_list <- lapply(ids, function(x)
        paste(x, formatC(seq_len(n_features),
                         digits = 2, flag = "0"), sep = "_"))

    out <- lapply(seq_len(n_types-1), FUN = function(x){
        randomLinkDF(l = id_list[-n_types][[x]],
                     r = id_list[-1L][[x]],
                     l_id = ids[-n_types][x],
                     r_id = ids[-1L][x],
                     p = (1-sparseness))})
    names(out) <- out_names
    asLinkMap(out)
}


#' Generate a random AnansiWeb, without any prior components
#' @description
#' called by randomWeb, not for user.
#' @importFrom Matrix rsparsematrix
#' @importFrom stats rnorm
#' @rdname randomAnansi
#' @noRd
#'
randomWebFull <- function(n_samp, n_x, n_y, density) {

    tableY <- matrix(data = rnorm(n_y * n_samp),
                     nrow = n_samp, ncol = n_y,
                     dimnames = list(
                         sample_id = paste0("sample_", seq_len(n_samp)),
                         y = paste0("y_", seq_len(n_y)))
    )
    tableX <- matrix(data = rnorm(n_x * n_samp),
                     nrow = n_samp, ncol = n_x,
                     dimnames = list(
                         sample_id = paste0("sample_", seq_len(n_samp)),
                         x = paste0("x_", seq_len(n_x)))
    )
    randomWebDic(tableY, tableX, density)
}

#' Generate a random AnansiWeb, only missing tables
#' @description
#' called by randomWeb, not for user.
#' @rdname randomAnansi
#' @noRd
#'
randomWebTab <- function(n_samp = 10, dictionary, metadata) {
    d <- dim(dictionary)
    tableY <- matrix(data = rnorm(d[1] * n_samp),
                     nrow = n_samp, ncol = d[1],
                     dimnames = list(
                         sample_id = paste0("sample_", seq_len(n_samp)),
                         y = rownames(dictionary))
    )
    names(dimnames(tableY))[2] <- names(dimnames(dictionary))[1]
    tableX <- matrix(data = rnorm(d[2] * n_samp),
                     nrow = n_samp, ncol = d[2],
                     dimnames = list(
                         sample_id = paste0("sample_", seq_len(n_samp)),
                         x = colnames(dictionary))
    )
    names(dimnames(tableX))[2] <- names(dimnames(dictionary))[2]
    metadata <- randomWebMetadata(tableY)
    # return AnansiWeb
    AnansiWeb( tableY = tableY, tableX = tableX,
               dictionary = dictionary, metadata = metadata)
}

#' Generate a random AnansiWeb, only missing dictionary.
#' @description
#' called by randomWeb, not for user.
#' @importFrom Matrix rsparsematrix
#' @rdname randomAnansi
#' @noRd
#'
randomWebDic <- function(tableY, tableX, density, metadata) {

    dictionary <- rsparsematrix(nrow = NCOL(tableY), ncol = NCOL(tableX),
                                density = density, rand.x = NULL,
                                dimnames = list(
                                    y = colnames(tableY),
                                    x = colnames(tableX)))
    names(dimnames(dictionary)) <- c(names(dimnames(tableY))[2L],
                                     names(dimnames(tableX))[2L])
    metadata <- randomWebMetadata(tableY)
    # return AnansiWeb
    AnansiWeb( tableY = tableY, tableX = tableX,
               dictionary = dictionary, metadata = metadata)
}

#' Generate random metadata for AnansiWeb
#' @description
#' called by randomWeb, not for user.
#' @param table a web table
#' @rdname randomAnansi
#' @noRd
#'
randomWebMetadata <- function(table){
    n_samples <- NROW(table)
    m <- DataFrame(
        cat_ab  = sample(c("a", "b"), n_samples, replace = TRUE),
        cat_XYZ = sample(c("X", "Y", "Z"), n_samples, replace = TRUE),
        num_norm  = rnorm(n_samples),
        num_unif  = runif(n_samples),
        row.names = paste0("sample_", seq_len(n_samples))
    )
    return(m)
}

#' Make a single df for a random AnansiLinkMap
#' @rdname randomAnansi
#' @description called by `randomLinkMap`, shouldn't be called by user.
#' @param l,r character vector of left, right features
#' @param l_id,r_id character scalar of left, right feature names
#' @param p proportion of connections to keep
#' @noRd
#'
randomLinkDF <- function(l, r, l_id, r_id, p){
    len <- length(l) * length(r)
    ind <- sort(sample(seq_len(len), size = ceiling(p*len)))
    out <- expand.grid(l, r, KEEP.OUT.ATTRS = FALSE)[ind,]
    names(out) <- c(l_id, r_id)
    out
}
