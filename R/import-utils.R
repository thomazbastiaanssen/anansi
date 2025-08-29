#' make a link data.frame for biobakery mapping files input
#' @param map `Character`, result from `readLines()` on uncompressed humann
#'     mapping files.
#' @returns a two-column `data.frame` that can be converted into an adjacency
#'     matrix used as input for `weaveWeb()`.
#' @seealso [weaveWeb()]
#' @export
#'
linkBiobakeryMap <- function(map){
    linkmap <- strsplit(x = map, split = "\t")
    linkmap_names <- vapply(linkmap, FUN = function(x) x[1], FUN.VALUE = "")
    linkmap <- lapply(linkmap, FUN = function(x) x[-1])

    linkmap <- data.frame(
        id.x = rep(linkmap_names, vapply(linkmap, length, 1)),
        id.y = unlist(linkmap, recursive = TRUE, use.names = FALSE)
    )
    return(linkmap)
}
