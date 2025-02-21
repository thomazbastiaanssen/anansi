#' Make a web
#' @description
#' Generate a biadjacency matrix linking the features between two tables.
#'
#' @param x name of feature type that will determine columns, corresponding to
#' `tableX`.
#' @param y name of feature type that will determine rows, corresponding to
#' `tableY`.
#' @param formula a formula of the form y ~ x, denoting desired output
#' format; assigns y to rows and columns to x. Equivalent to using \code{x} and
#' \code{y} arguments.
#' @param link a data.frame with two columns, as named in \code{x} and \code{y}.
#' Sharing a row indicates a link between features. See \code{ec2ko} and
#' \code{ec2cpd} for examples. Optionally a list with two such data.frames.
#' @param x.ids a character vector of features of type \code{x}, which should be
#' included. If left NULL (default), keep all features.
#' @param y.ids a character vector of features of type \code{y}, which should be
#' included. If left NULL (default), keep all features.
#' @param ... further arguments to be passed to or from methods. Not used.
#'
#' @returns a sparse binary biadjacency matrix, with features from \code{y} as
#' rows and features from \code{x} as columns.
#'
#' @export
#' @examples
#' # Basic usage
#' web(x = "ko", y = "ec", link = ec2ko)
#' web(ec ~ cpd, link = ec2cpd)
#'
#' # A wrapper is available for kegg ko, ec and cpd data
#' generic      <- web(cpd ~ ko, link = list(ec2ko, ec2cpd))
#' kegg_wrapper <- kegg_web(cpd ~ ko)
#'
#' identical(generic, kegg_wrapper)
#'
#' # The following are equivalent:
#' a <- web(ko ~ cpd, link = list(ec2ko, ec2cpd))
#' b <- web(cpd ~ ko, link = list(ec2ko, ec2cpd))
#'
#' identical(a, Matrix::t(b))
#'
web <- function(x, ...) UseMethod("web")

#' @rdname web
#' @export
#'
kegg_web <- function(x, ...){
  web(x, link = list(anansi::ec2ko, anansi::ec2cpd), ...)
}

#' @rdname web
#' @export
#'
web.default <- function(x, y, link, x.ids = NULL, y.ids = NULL, ...){
  terms <- c(x, y)
  stopifnot("both 'x' and 'y' terms must be provided as character" =
              is(terms, "character") && length(terms) == 2L)

  if(!is.null(dim(link)) ) df.names <- colnames(link)

  if( is.null(dim(link)) ) {
    df.names <- unique(unlist(lapply(link, colnames)))
    cn.1 <- colnames(link[[1L]]); cn.2 <- colnames(link[[2L]])
    if(all(cn.1 %in% terms)) link <- link[[1L]]
    if(all(cn.2 %in% terms)) link <- link[[2L]]
    }

  stopifnot("both 'x' and 'y' terms must be found as colnames in 'link'" =
              all(terms %in% df.names) )

  if(!is.null(dim(link)) ) {

    link <- link[,terms]

    if(!is.null(x.ids))  link <- link[link[,1L] %in% x.ids,]
    if(!is.null(y.ids))  link <- link[link[,2L] %in% y.ids,]

    df_to_sparse_biadjacency_matrix(link)

  } else if (is.null(dim(link))) {
    i <- intersect(cn.1, cn.2)
    stopifnot("data.frames in 'link' must share a colname" = length(i) == 1L)

    if(all(cn.1 %in% c(i, terms[1]))){

      x.df <- link[[1L]][, c(i, terms[1L])]
      y.df <- link[[2L]][, c(i, terms[2L])]

    } else if(all(cn.1 %in% c(i, terms[2]))){

      x.df <- link[[2L]][, c(i, terms[1L])]
      y.df <- link[[1L]][, c(i, terms[2L])]
    }
      web_from_2_dfs(x.df, y.df, x.ids, y.ids)
    }
  }

#' @rdname web
#' @export
#'
web.formula <- function(formula, link, x.ids = NULL, y.ids = NULL, ...){

  if (missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")

  terms <- all.vars(formula)

  if(!is.null(dim(link)) ) df.names <- colnames(link)
  if( is.null(dim(link)) ) df.names <- unique(unlist(lapply(link, colnames)))

  if (sum(terms %in% df.names) != 2L)
    stop("Exactly both colnames should be used in 'formula'")

  web.default(x = terms[2], y = terms[1], link, x.ids, y.ids)
}


#' @importFrom Matrix tcrossprod
#'
web_from_2_dfs <- function(x.df, y.df, x.ids, y.ids){
  if(!is.null(x.ids)) { x.df <- x.df[x.df[,2L] %in% x.ids,] }
  if(!is.null(y.ids)) { y.df <- y.df[y.df[,2L] %in% y.ids,] }

  i <- intersect(x.df[,1L], y.df[,1L])
  x.df <- x.df[x.df[,1L] %in% i,]
  y.df <- y.df[y.df[,1L] %in% i,]
  x.mat <- df_to_sparse_biadjacency_matrix(x.df)
  y.mat <- df_to_sparse_biadjacency_matrix(y.df)

  tcrossprod(y.mat, x.mat, boolArith = TRUE)
}

#' @importFrom igraph vertex_attr<- vertex_attr graph_from_data_frame as_biadjacency_matrix
#'
df_to_sparse_biadjacency_matrix <- function(x){
  x.g <- graph_from_data_frame(x, directed = FALSE)
  vertex_attr(x.g, name = "type") <- vertex_attr(x.g, "name") %in% x[,1L]
  as_biadjacency_matrix(x.g, sparse = TRUE)
}

