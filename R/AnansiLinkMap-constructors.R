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
        paste(x, formatC(seq_len(n_features), digits = 2, flag = "0"), sep = "_"))

    out <- lapply(seq_len(n_types-1), FUN = function(x){
        randomLinkDF(l = id_list[-n_types][[x]],
                     r = id_list[-1L][[x]],
                     l_id = ids[-n_types][x],
                     r_id = ids[-1L][x],
                     p = (1-sparseness))})
    names(out) <- out_names
    asLinkMap(out)
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
