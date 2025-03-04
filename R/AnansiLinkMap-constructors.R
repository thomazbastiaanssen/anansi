#' Make and anansiLinkMap
#' @name LinkMap
#' @rdname LinkMap
#' @aliases anansiLinkMap asLinkMap
#' @description
#' \code{asLinkMap()} constructs an \code{anansiLinkMap} object from a validly
#' formed data frame or list of such data.frames.
#' @param x \code{any} object, most likely \code{list} of data frames.
#' @export
#' @seealso \itemize{
#' \item \code{\link{kegg_link}}: for an example of valid input.
#' \item \code{\link{anansiLinkMap-class}}: for class.
#' \item \code{\link{anansiLinkMap-methods}} for methods.
#'}
#' @examples
#' asLinkMap( kegg_link() )
#'
asLinkMap <- function(x) {
  if(validLinkDF(x)) x <- list(link = x)

  linkMap <- new("anansiLinkMap", x)
  validObject(linkMap)

  return(linkMap)
}
