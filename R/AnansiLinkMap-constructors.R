#' Make an AnansiLinkMap
#' @name AnansiLinkMap
#' @rdname AnansiLinkMap
#' @aliases LinkMap asLinkMap
#' @description
#' Construct an `AnansiLinkMap` object from a validly shaped data frame or
#' list of such data.frames.
#' @param x `any` object, most likely `list` of data frames.
#' @export
#' @seealso \itemize{
#' \item [kegg_link()]: for an example of valid input.
#' \item [AnansiLinkMap-class()]: for class.
#' \item [AnansiLinkMap-methods()] for methods.
#'}
#' @examples
#' AnansiLinkMap( kegg_link( ) )
#'
AnansiLinkMap <- function(x) {
  if(validLinkDF(x)) x <- list(link = x)

  linkMap <- new("AnansiLinkMap", x)
  validObject(linkMap)

  return(linkMap)
}

#' @rdname AnansiLinkMap
#' @export
#'
asLinkMap <- AnansiLinkMap

#' @rdname AnansiLinkMap
#' @description
#' Helper function that takes an AnansiLinkMap and returns a sparse biadjacency
#' Matrix with link df names as rownames and id names as colnames. Called
#' internally.
#' @returns a sparse biadjacency Matrix with link df names as rownames and id
#'     names as colnames
#' @importFrom Matrix sparseMatrix
#'
linkMatrix <- function(x){
    i <- factor(rep(rownames(x), each = 2))
    j <- factor(unlist(names(x), use.names = FALSE))
    sparseMatrix(i, j, dimnames = dimnames(x))
}

