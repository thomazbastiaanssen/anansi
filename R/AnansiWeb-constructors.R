#' @name AnansiWeb
#' @rdname AnansiWeb
#' @param tableY,tableX A table containing features of interest. Rows should be
#'     samples and columns should be features. Y and X refer to the position of
#'     the features in a formula: Y ~ X.
#' @param dictionary A binary adjacency matrix of class `Matrix`, or coercible
#'     to `Matrix`
#' @param ... further arguments.
#' @seealso \itemize{
#'  \item [kegg_link()]: For examples of input for link argument.
#' }
#' @returns an `AnansiWeb` object, with sparse binary biadjacency matrix
#' with features from `y` as rows and features from `x` as columns in
#' `dictionary` slot.
#' @usage
#' ## Constructor for `AnansiWeb` objects
#' AnansiWeb(tableX, tableY, dictionary, metadata = list(), ...)
#' @examples
#'
#' # Use AnansiWeb() to consrtuct an AnansiWeb object from components:
#' tX <- `dimnames<-`(replicate(5, (rnorm(36))),
#'     value = list(
#'         as.character(seq_len(36)),
#'         letters[1:5]
#'     )
#' )
#' tY <- `dimnames<-`(replicate(3, (rnorm(36))),
#'     value = list(
#'         as.character(seq_len(36)),
#'         LETTERS[1:3]
#'     )
#' )
#'
#' d <- matrix(TRUE,
#'     nrow = NCOL(tY), ncol = NCOL(tX),
#'
#'     # Note: Dictionary should have named dimensions
#'     dimnames = list(
#'         y_names = colnames(tY),
#'         x_names = colnames(tX)
#'     )
#' )
#' web <- AnansiWeb(tableX = tX, tableY = tY, dictionary = d)
#' @param metadata `list` of metadata. Optional.
#' @importFrom Matrix Matrix drop0
#' @importFrom S4Vectors DataFrame
#' @export
#'
AnansiWeb <- function(tableX, tableY, dictionary, metadata = list(), ...) {
    # coerce
    if (!is(dictionary, "Matrix")) {
        dictionary <-
            drop0(Matrix(dictionary, sparse = TRUE))
    }
    if (!is(tableX, "matrix")) tableX <- as.matrix(tableX)
    if (!is(tableY, "matrix")) tableY <- as.matrix(tableY)

    # check validity
    stopifnot(
        "'tableX' and 'tableY' need same number of rows (observations)" = NROW(
            tableX
        ) ==
            NROW(tableY)
    )
    stopifnot(
        "cols in 'tableY' need same amount as rows in dictionary" = NCOL(
            tableY
        ) ==
            NROW(dictionary)
    )
    stopifnot(
        "cols in 'tableX' need same amount as rows in dictionary" = NCOL(
            tableX
        ) ==
            NCOL(dictionary)
    )
    if (
        is.null(names(dimnames(dictionary))) ||
            any(names(dimnames(dictionary)) %in% "")
    ) {
        warning("Dimnames of 'dictionary' were missing; Assigned 'y' and 'x'.")
        names(dimnames(dictionary)) <- c("y", "x")
    }

    if (!inherits(metadata, "list")) {
        metadata <-
            list(metadata = as.data.frame(metadata))
    }
    # return AnansiWeb
    new(
        "AnansiWeb",
        tableY = tableY,
        tableX = tableX,
        dictionary = dictionary,
        metadata = metadata
    )
}
