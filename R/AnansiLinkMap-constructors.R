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

#' Make a random anansiLinkMap
#' @rdname LinkMap
#' @export
#'
randomLinkMap <- function(n_ids = 6, n_feat = 100, p = 1){
    n_ids <- max(min(n_ids, 26), 2)
    ids <- letters[seq_len(n_ids)]
    out_names <- paste0(ids[-n_ids], "2",ids[-1L])
    id_list <- lapply(ids, function(x)
        paste(x, formatC(seq_len(n_feat), digits = 2, flag = "0"), sep = "_"))

    out <- lapply(seq_len(n_ids-1), FUN = function(x){
        randomLinkDF(l = id_list[-n_ids][[x]],
                     r = id_list[-1L][[x]],
                     l_id = ids[-n_ids][x],
                     r_id = ids[-1L][x],
                     p = p)})
    names(out) <- out_names
    asLinkMap(out)
}

#' Make a single df for a random anansiLinkMap
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
