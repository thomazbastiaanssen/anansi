#' Coercion methods for AnansiWeb
#' @description
#' Coerce AnansiWeb to and from other object types.
#' @name coerceAnansi
#' @examples
#' # Create a random web
#' web <- randomWeb()
#'
#' # To list
#' as.list(web)
#'
#' # To MultiAssayExperiment
#' asMAE(web)
#'
NULL

#' @rdname coerceAnansi
#' @aliases as.list.AnansiWeb
#' @inheritParams base::as.list
#' @export
#'
setMethod("as.list", c(x = "AnansiWeb"), function(x, ...) as(x, "list"))

#' @rdname coerceAnansi
#' @aliases as.MAE as.MultiAssayExperiment asMultiAssayExperiment
#' @export
#'
asMAE <- function(x) as(x, "MultiAssayExperiment")

#' @importFrom methods as
#' @export
#'
setAs(from = "AnansiWeb", to = "list", def = function(from) {
    out <- c(list(
        tableY = from@tableY, tableX = from@tableX,
        dictionary = from@dictionary
    ), from@metadata)
    names(out)[c(1L, 2L)] <- names(from)
    out
})

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#'
setAs(from = "AnansiWeb", to = "MultiAssayExperiment", def = function(from) {
    tY <- t(tableY(from))
    tX <- t(tableX(from))
    to_exp <- ExperimentList(
        y = SummarizedExperiment(tY),
        x = SummarizedExperiment(tX)
    )
    names(to_exp) <- names(from)

    to_md <- list(dictionary = from@dictionary)
    to_cd <- metadata(from, simplify = TRUE)

    MultiAssayExperiment(
        experiments = to_exp,
        metadata = to_md,
        colData = DataFrame(to_cd)
    )
})

#' @description Convert MultiFactor to list
#' @rdname coerceAnansi
#' @aliases as.list.MultiFactor
#' @inheritParams BiocGenerics::as.list
#' @returns a named list of character vectors (Default) or integers
#' (`use.names = FALSE`).
#' @export
#'
setMethod("as.list", c(x = "MultiFactor"), function(x, ..., use.names = TRUE) {
    as.list.MultiFactor(x, ..., use.names)
})

#' @export
#' @rdname coerceAnansi
#' @param use.names `Logical scalar`, whether output list should contain
#'     character (Default) or integer data frame. If `FALSE`, returns
#'     `unfactor(x)`.
#' @seealso [unfactor()]
#' @examples
#' x <- as.list(randomMultiFactor())
#' identical(x, as.list(MultiFactor(x)))
#'
as.list.MultiFactor <- function(x, use.names = TRUE) {
    ifelse(
        use.names,
        yes = return(unfactor(x)),
        no  = return(`names<-`(x@index, rownames(x)))
    )
}
