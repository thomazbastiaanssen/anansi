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
#' @param tableY A table containing features of interest. Rows should be samples
#' and columns should be features. The Y and X refer to the position of the
#' features in a formula: Y ~ X.
#' @param tableX A table containing features of interest. Rows should be samples
#' and columns should be features. The Y and X refer to the position of the
#' features in a formula: Y ~ X.
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
#' kegg_wrapper <- kegg_web( cpd ~ ko )
#'
#' identical(generic, kegg_wrapper)
#'
#' # The following are equivalent to transposition:
#' a <- anansi:::get_dict( web(ko ~ cpd, link = list(ec2ko, ec2cpd)) )
#' b <- anansi:::get_dict( web(cpd ~ ko, link = list(ec2ko, ec2cpd)) )
#'
#' identical(a, t(b))
#'
web <- function(x, ...) UseMethod("web")

#' @rdname web
#' @export
#'
kegg_web <- function(x, ...){
  web(x, link = list(anansi::ec2ko, anansi::ec2cpd), ...)
}

#' @rdname web
#' @importFrom Matrix Matrix
#' @export
#'
web.default <- function(x, y, link = NULL, tableX = NULL, tableY = NULL, ...){
  terms <- c(x, y)
  stopifnot("both 'x' and 'y' terms must be provided as character" =
              is(terms, "character") && length(terms) == 2L)

  if(identical(link,"none")) return(web_missing_link(tableX, tableY))
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
  if(is.null(tableX) && is.null(tableY)) return(
    new("anansiWeb",
      tableY     = matrix(ncol = NROW(d), dimnames = list(NULL, rownames(d))),
      tableX     = matrix(ncol = NCOL(d), dimnames = list(NULL, colnames(d))),
      dictionary = d)) else return(
        new("anansiWeb",
        tableY     = tableY[,rownames(d)],
        tableX     = tableX[,colnames(d)],
        dictionary = d))
}

#' @rdname web
#' @export
#'
web.formula <- function(formula, link = NULL, tableX = NULL, tableY = NULL,
                        ...){

  if (missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")

  terms <- all.vars(formula)
  if(is.null(link)) return(
    web.default(x = terms[2], y = terms[1], link, tableX, tableY)
    )

  if(!is.null(dim(link)) ) df.names <- colnames(link)
  if( is.null(dim(link)) ) df.names <- unique(unlist(lapply(link, colnames)))

  if (sum(terms %in% df.names) != 2L)
    stop("Exactly both colnames should be used in 'formula'")

  web.default(x = terms[2], y = terms[1], link, tableX, tableY)
}

#' Produce a biadjacency matrix given tables and a dictionary
#' @description calculates a biadjacency matrix for the cases where
#' \code{link} is a single - or a \code{list} of two - \code{data.frame}(s).
#' @param link a \code{data.frame} or \code{list} of two compatible ones.
#' @param terms a length 2 character vector, naming the x and y terms in order
#' @param tableY A table containing features of interest. Rows should be samples
#' and columns should be features. The Y and X refer to the position of the
#' features in a formula: Y ~ X.
#' @param tableX A table containing features of interest. Rows should be samples
#' and columns should be features. The Y and X refer to the position of the
#' features in a formula: Y ~ X.
#' @param link_is_list a boolean, is a \code{list} (or \code{FALSE}: a
#' \code{data.frame})
#' @returns a sparse boolean biadjacency matrix, for use in main anansi workflow
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
    i <- intersect(cn.1, cn.2)
    stopifnot("data.frames in 'link' must share a colname" = length(i) == 1L)

    if(all(cn.1 %in% c(i, terms[1]))){

      x.df <- link[[1L]][, c(i, terms[1L])]
      y.df <- link[[2L]][, c(i, terms[2L])]

    } else if(all(cn.1 %in% c(i, terms[2]))){

      x.df <- link[[2L]][, c(i, terms[1L])]
      y.df <- link[[1L]][, c(i, terms[2L])]

    }
    web_from_2_dfs(x.df, y.df, colnames(tableX), colnames(tableY))
  }
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

#' Make a full web; for all vs all association testing
#' @description
#' Make a fully TRUE biadjacency matrix with dimensions of the two input tables.
#' @param tableX a matrix of features of table \code{X}.
#' @param tableY a matrix of features of table \code{Y}.
#' @returns An \code{anansiWeb} object with both tables and a fully \code{TRUE}
#' (non-sparse) matrix from the \code{Matrix} package.
#' @importFrom Matrix Matrix
#'
web_missing_link <- function(tableX, tableY) {

  d <- Matrix(
    data = TRUE,
    nrow = NROW(tableY),
    ncol = NROW(tableX),
    dimnames = list(sort(colnames(tableY)), sort(colnames(tableX)))
    )

  new("anansiWeb",
      tableY     = tableY[,rownames(d)],
      tableX     = tableX[,colnames(d)],
      dictionary = d)

}

