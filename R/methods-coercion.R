#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases as.list.AnansiWeb coerce,AnansiWeb-list
#' @usage NULL
#' @export
#'
method(as.list, AnansiWeb)  <- function(x) {
    out <- c(
        list(
            tableY = x@tableY,
            tableX = x@tableX,
            dictionary = x@dictionary
        ),
        x@metadata
    )
    names(out)[c(1L, 2L)] <- names(x)
    out
}

#' @name AnansiWeb
#' @rdname AnansiWeb
#' @param row.names,optional Ignored, for S4 generic. See ?base::as.data.frame.
#' @aliases as.data.frame.AnansiWeb-method coerce,AnansiWeb-data.frame
#' @usage NULL
#' @export
#'
method(as.data.frame, AnansiWeb) <- function(x) {
    cbind(
        tableY(x),
        tableX(x),
        metadata(x, simplify = TRUE)
    )
}

#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases as.MAE as.MultiAssayExperiment asMultiAssayExperiment
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#'
asMAE <- function(x)  {
    tY <- t(tableY(x))
    tX <- t(tableX(x))
    to_exp <- ExperimentList(
        y = SummarizedExperiment(tY),
        x = SummarizedExperiment(tX)
    )
    names(to_exp) <- names(x)

    to_md <- list(dictionary = dictionary(x))
    to_cd <- metadata(x, simplify = TRUE)

    MultiAssayExperiment::MultiAssayExperiment(
        experiments = to_exp,
        metadata = to_md,
        colData = DataFrame(to_cd)
    )
}

#' @name MultiFactor
#' @rdname MultiFactor
#' @aliases as.list.MultiFactor
#' @returns a named list of character vectors (Default) or integers
#' (`use.names = FALSE`).
#' @param use.names `Logical scalar`, whether output list should contain
#'     character (Default) or integer data frame. If `FALSE`, returns
#'     `unfactor(x)`.
#' @seealso [unfactor()]
#' @export
#'
method(as.list, MultiFactor) <- function(x, use.names = TRUE) {
    ifelse(
        use.names,
        yes = return(unfactor(x)),
        no = return(`names<-`(x@index, rownames(x)))
    )}
