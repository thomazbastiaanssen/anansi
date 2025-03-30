#' @rdname AnansiWeb
#' @aliases as.list.AnansiWeb
#' @export
#'
setMethod("as.list", c(x = "AnansiWeb"), function(x, ...) as(x, "list"))

#' @rdname AnansiWeb
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

    to_md <- list(dictionary = dictionary(from))
    to_cd <- metadata(from, simplify = TRUE)

    MultiAssayExperiment::MultiAssayExperiment(
        experiments = to_exp,
        metadata = to_md,
        colData = DataFrame(to_cd)
    )
})

#' @rdname MultiFactor
#' @aliases as.list.MultiFactor
#' @returns a named list of character vectors (Default) or integers
#' (`use.names = FALSE`).
#' @export
#'
setMethod("as.list", c(x = "MultiFactor"), function(x, ..., use.names = TRUE) {
    as.list.MultiFactor(x, ..., use.names)
})

#' @export
#' @method as.list MultiFactor
#' @rdname MultiFactor
#' @param use.names `Logical scalar`, whether output list should contain
#'     character (Default) or integer data frame. If `FALSE`, returns
#'     `unfactor(x)`.
#' @seealso [unfactor()]
#'
as.list.MultiFactor <- function(x, ..., use.names = TRUE) {
    ifelse(
        use.names,
        yes = return(unfactor(x)),
        no  = return(`names<-`(x@index, rownames(x)))
    )
}
