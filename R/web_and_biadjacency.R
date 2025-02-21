#'
#' web(x = "ko", y = "cpd", x.df = ec2ko, y.df = ec2cpd)
#' web(ec ~ cpd, x.df = ec2cpd)
#' web(ko ~ cpd, y.df = ec2cpd, x.df = ec2ko)
#' kegg_web(ec~ko)
#'
web <- function(x, ...) UseMethod("web")

#' @export
#'
kegg_web <- function(x, ...){
  web(x, x.df = anansi::ec2ko, y.df = anansi::ec2cpd, ...)
}

#' @exportS3Method
#'
web.default <- function(x, y, x.df, y.df = NULL, x.ids = NULL, y.ids = NULL){
  terms <- c(x, y)
  stopifnot("both 'x' and 'y' terms must be provided as character" =
              is(terms, "character") & length(terms) == 2L)

  stopifnot("both 'x' and 'y' terms must be found as colnames for x.df; y.df" =
              all(terms %in% c(colnames(x.df), colnames(y.df))))

  if( is.null(y.df) | all(terms %in% colnames(x.df)) ){
    x.df <- x.df[,terms]

    if(!is.null(x.ids))  x.df <- x.df[x.df[,1L] %in% x.ids,]
    if(!is.null(y.ids))  x.df <- y.df[y.df[,2L] %in% y.ids,]

      df_to_sparse_biadjacency_matrix(x.df)

  } else if( all(terms %in% colnames(y.df)) ){

    y.df <- y.df[,terms]

    if(!is.null(x.ids))  y.df <- y.df[y.df[,1L] %in% x.ids,]
    if(!is.null(y.ids))  y.df <- y.df[y.df[,2L] %in% y.ids,]

    df_to_sparse_biadjacency_matrix(y.df)

  } else {
    i <- intersect(colnames(x.df), colnames(y.df))
    stopifnot("x.df and y.df must share a colname" = length(i) == 1L)
    if(terms[1L] %in% colnames(x.df)){
      web_from_2_dfs(x.df[, c(i, terms[1L])],
                     y.df[, c(i, terms[2L])],
                     x.ids, y.ids)

    } else if(terms[1L] %in% colnames(y.df)){
      web_from_2_dfs(y.df[, c(i, terms[1L])],
                     x.df[, c(i, terms[2L])],
                     x.ids, y.ids)
    }
  }
  }

#' @exportS3Method
#'
web.formula <- function(formula, x.df, y.df = NULL, y.ids = NULL, x.ids = NULL){

  if (missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")

  terms <- all.vars(formula)
  df.names <- c(colnames(x.df), colnames(y.df))

  if (sum(terms %in% df.names) != 2L)
    stop("Exactly both colnames should be used in 'formula'")

  web.default(x = terms[2], y = terms[1], x.df, y.df, x.ids, y.ids)
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

