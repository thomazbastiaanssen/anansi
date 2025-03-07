#' Make an AnansiWeb
#' @name AnansiWeb
#' @rdname AnansiWeb
#' @description
#' Generate a biadjacency matrix, linking the features between two tables.
#' Return an `AnansiWeb` object which contains all three.
#'
#' `weaveWeb()` is for general use and has flexible default settings.
#'
#' `weaveKEGG()` is a wrapper that sets `link` to `kegg_link()`.
#' All variants are special cases of `weaveWeb()`.
#'
#' `AnansiWeb()` constructs an `AnansiWeb` object from three tables.
#'
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
#'  \item [AnansiWeb-methods()]: For utility functions to get and set.
#'  \item [kegg_link()]: For examples of input for link argument.
#'  \item [getWeb()]: For
#'  [MultiAssayExperiment::MultiAssayExperiment()] methods.
#' }
#'
#' @returns an `AnansiWeb` object, with sparse binary biadjacency matrix
#' with features from `y` as rows and features from `x` as columns in
#' `dictionary` slot.
#' @examples
#' # Basic usage
#' weaveWeb(cpd ~ ko, link = kegg_link())
#' weaveWeb(x = "ko", y = "ec", link = ec2ko)
#' weaveWeb(ec ~ cpd, link = ec2cpd)
#'
#' # use AnansiWeb() to constuct an AnansiWeb object from components:
#'
#' tX <- `colnames<-`(replicate(5, (rnorm(36))), letters[1:5])
#' tY <- `colnames<-`(replicate(3, (rnorm(36))), LETTERS[1:3])
#' d <- matrix(TRUE, nrow = NCOL(tY), ncol = NCOL(tX),
#'             dimnames = list(y = colnames(tY), x = colnames(tX)))
#'
#' AnansiWeb(tableX = tX, tableY = tY, dictionary = d)
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
NULL

#' @rdname AnansiWeb
#' @order 0
#' @usage NULL
#' @export
#'
weaveWeb <- function(x, ...) UseMethod("weaveWeb")

#' @rdname AnansiWeb
#' @importFrom Matrix Matrix
#' @order 2
#' @export
#'
weaveWeb.default <- function(x, y, link = NULL, tableX = NULL, tableY = NULL, ...){
  terms <- c(x, y)
  stopifnot("both 'x' and 'y' terms must be provided as character" =
              is(terms, "character") && length(terms) == 2L)

  if(identical(link,"none")) return(web_missing_link(tableX, tableY, x, y))

  link_is_list <- is.list(link) && !is.data.frame(link)

  if(link_is_list) {
    df.names <- unique(unlist(lapply(link, colnames)))

    stopifnot("both 'x' and 'y' terms must be found as colnames in 'link'" =
                all(terms %in% df.names) )

    cn.1 <- colnames(link[[1L]]); cn.2 <- colnames(link[[2L]])
    if(all(cn.1 %in% terms)) link <- link[[1L]]
    if(all(cn.2 %in% terms)) link <- link[[2L]]
    link_is_list <- is.list(link) && !is.data.frame(link)
  }
  d <- deliver_web_cases(link, terms, tableX, tableY, link_is_list)
  names(dimnames(d)) <- c(y, x)
  if(is.null(tableX) && is.null(tableY))
    return(
      AnansiWeb(
          tableY     = matrix(ncol = NROW(d),
                              dimnames = list(NULL, rownames(d))),
          tableX     = matrix(ncol = NCOL(d),
                              dimnames = list(NULL, colnames(d))),
          dictionary = d,
          metadata   = DataFrame(metadata))
    ) else
      return(
        AnansiWeb(
            tableY     = as.matrix(tableY)[,rownames(d)],
            tableX     = as.matrix(tableX)[,colnames(d)],
            dictionary = d,
            metadata   = DataFrame(metadata))
        )
}

#' @rdname AnansiWeb
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

  if(!is.null(dim(link)) ) df.names <- colnames(link)
  if( is.null(dim(link)) ) df.names <- unique(unlist(lapply(link, colnames)))

  if (sum(terms %in% df.names) != 2L)
    stop("Exactly both colnames should be used in 'formula'")

  weaveWeb.default(x = terms[2], y = terms[1], link, tableX, tableY)
}

#' @rdname AnansiWeb
#' @param dictionary A binary adjacency matrix of class `Matrix`, or
#' coercible to `Matrix`
#' @param metadata `DataFrame`, or coerible to `DataFrame`. Optional.
#' @importFrom Matrix Matrix drop0
#' @importFrom S4Vectors DataFrame
#' @export
#'
AnansiWeb <- function(tableX, tableY, dictionary, metadata = NULL, ...) {
  # coerce
  if(!is(dictionary, "Matrix")) dictionary <-
      drop0(Matrix(dictionary, sparse = TRUE))
  if(!is(tableX, "matrix")) tableX <- as.matrix(tableX)
  if(!is(tableY, "matrix")) tableY <- as.matrix(tableY)

  # check validity
  stopifnot("'tableX' and 'tableY' need same number of rows (observations)" =
            NROW(tableX) == NROW(tableY))
  stopifnot("cols in 'tableY' need same amount as rows in dictionary" =
              NCOL(tableY) == NROW(dictionary))
  stopifnot("cols in 'tableX' need same amount as rows in dictionary" =
              NCOL(tableX) == NCOL(dictionary))
  if( is.null( names(dimnames(dictionary)) ) ||
      any( names(dimnames(dictionary)) %in% "")) {
    warning("Dimnames of 'dictionary' were missing; Assigned 'y' and 'x'.")
       names(dimnames(dictionary)) <- c("y", "x")
    }

  # return AnansiWeb
  new("AnansiWeb",
        tableY     = tableY,
        tableX     = tableX,
        dictionary = dictionary,
        metadata   = DataFrame(metadata))
    }

#' @rdname AnansiWeb
#' @export
#'
weaveKEGG <- function(x, ...) weaveWeb(x, link = kegg_link(), ...)

#' Generate a random AnansiWeb or AnansiLinkMap
#' @rdname randomAnansi
#' @param tableY,tableX A table containing features of interest. Rows should be
#'     samples and columns should be features. Y and X refer to the position of
#'     the features in a formula: Y ~ X.
#' @param dictionary A binary adjacency matrix of class `Matrix`, or
#' coercible to `Matrix`
#' @export
#' @seealso [AnansiWeb()], [AnansiLinkMap()]
#' @examples
#' # Make a random AnansiWeb object
#' randomWeb()
#' randomLinkMap()
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

###############################################################################
###############################################################################

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
#' @param link_is_list a boolean, is a `list` (or `FALSE`: a
#' `data.frame`)
#' @returns a sparse boolean biadjacency matrix, for use in main anansi workflow
#' @noRd
#'
deliver_web_cases <- function(link, terms, tableX, tableY, link_is_list){
  if(is.data.frame(link)) {
    df.names <- colnames(link)

    stopifnot("both 'x' and 'y' terms must be found as colnames in 'link'" =
                all(terms %in% df.names) )

    link <- link[,terms]

    if(!is.null(tableX))  link <- link[link[,1L] %in% colnames(tableX),]
    if(!is.null(tableY))  link <- link[link[,2L] %in% colnames(tableY),]

    return( df_to_sparse_biadjacency_matrix(link) )

  } else if (link_is_list) {
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

    web_from_2_dfs(x.df, y.df, colnames(tableX), colnames(tableY))
  }
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
#' @param x,y `character scalar` names of x & y terms
#' @returns
#' An `AnansiWeb` object with both tables and a fully `TRUE`
#' (non-sparse) matrix from the `Matrix` package.
#' @importFrom Matrix Matrix
#' @noRd
#'
web_missing_link <- function(tableX, tableY, x, y, metadata = NULL) {

  d <- Matrix(
    data = TRUE,
    nrow = NCOL(tableY),
    ncol = NCOL(tableX),
    dimnames = list(sort(colnames(tableY)), sort(colnames(tableX)))
    )
  names(dimnames(d)) <- c(y, x)

  AnansiWeb(
      tableY     = as.matrix(tableY)[,rownames(d)],
      tableX     = as.matrix(tableX)[,colnames(d)],
      dictionary = d,
      metadata = DataFrame(metadata))

}
