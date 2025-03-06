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
#' @examples asLinkMap( kegg_link() )
#' @export
asLinkMap <- AnansiLinkMap
