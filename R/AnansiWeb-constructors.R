#' Make an AnansiWeb
#' @name AnansiWeb
#' @rdname AnansiWeb
#' @description
#' `AnansiWeb()` constructs an `AnansiWeb` object from three tables.
#'
#' @param link One of the following:
#' \itemize{
#'  \item `Character scalar` with value `"none"`.
#'  \item `data.frame` with two columns
#'  \item `list` with two such `data.frame`s.
#' }
#' @param tableY,tableX A table containing features of interest. Rows should be
#'     samples and columns should be features. Y and X refer to the position of
#'     the features in a formula: Y ~ X.
#' @param ... further arguments.
#' @details
#' If the `link` argument is `"none"`, all features will be considered
#' linked. If one or more `data.frame`s, colnames should be as specified in
#' `x` and `y`.
#' @seealso \itemize{
#'  \item [AnansiWeb-methods()]: For utility functions to get and set.
#'  \item [kegg_link()]: For examples of input for link argument.
#'  \item [getWeb()]: For
#'  [MultiAssayExperiment::MultiAssayExperiment()] methods.
#' }
#'
#' @returns an `AnansiWeb` object, with sparse binary biadjacency matrix
#' with features from `y` as rows and features from `x` as columns in
#' `dictionary` slot.
#' @examples
#'
#' # use AnansiWeb() to constuct an AnansiWeb object from components:
#'
#' tX <- `colnames<-`(replicate(5, (rnorm(36))), letters[1:5])
#' tY <- `colnames<-`(replicate(3, (rnorm(36))), LETTERS[1:3])
#' d <- matrix(TRUE, nrow = NCOL(tY), ncol = NCOL(tX),
#'             dimnames = list(y = colnames(tY), x = colnames(tX)))
#'
#' AnansiWeb(tableX = tX, tableY = tY, dictionary = d)
#'
NULL

#' @rdname AnansiWeb
#' @param dictionary A binary adjacency matrix of class `Matrix`, or
#' coercible to `Matrix`
#' @param metadata `list` of metadata. Optional.
#' @importFrom Matrix Matrix drop0
#' @importFrom S4Vectors DataFrame
#' @export
#'
AnansiWeb <- function(tableX, tableY, dictionary, metadata = list(), ...) {
  # coerce
  if(!is(dictionary, "Matrix")) dictionary <-
      drop0(Matrix(dictionary, sparse = TRUE))
  if(!is(tableX, "matrix")) tableX <- as.matrix(tableX)
  if(!is(tableY, "matrix")) tableY <- as.matrix(tableY)

  # check validity
  stopifnot("'tableX' and 'tableY' need same number of rows (observations)" =
            NROW(tableX) == NROW(tableY))
  stopifnot("cols in 'tableY' need same amount as rows in dictionary" =
              NCOL(tableY) == NROW(dictionary))
  stopifnot("cols in 'tableX' need same amount as rows in dictionary" =
              NCOL(tableX) == NCOL(dictionary))
  if( is.null( names(dimnames(dictionary)) ) ||
      any( names(dimnames(dictionary)) %in% "")) {
    warning("Dimnames of 'dictionary' were missing; Assigned 'y' and 'x'.")
       names(dimnames(dictionary)) <- c("y", "x")
    }

  if(!inherits(metadata, "list")) metadata <-
    list(metadata = as.data.frame(metadata))
  # return AnansiWeb
  new("AnansiWeb",
        tableY     = tableY,
        tableX     = tableX,
        dictionary = dictionary,
        metadata   = metadata)
    }
