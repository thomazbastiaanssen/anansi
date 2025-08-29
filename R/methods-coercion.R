#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases as.list.AnansiWeb coerce,AnansiWeb-list
#' @usage NULL
#' @export
#'
S7::method(as.list, AnansiWeb)  <- function(x) {
    out <- S7::props(x)
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
S7::method(as.data.frame, AnansiWeb) <- function(x) {
    cbind(
        x@tableY,
        x@tableX,
        x@metadata
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

    tY <- t(x@tableY)
    tX <- t(x@tableX)
    to_exp <- ExperimentList(
        y = SummarizedExperiment::SummarizedExperiment(tY),
        x = SummarizedExperiment::SummarizedExperiment(tX)
    )
    names(to_exp) <- names(x)

    to_md <- list(dictionary = x@dictionary)
    to_cd <- x@metadata

    if(prod(dim.data.frame(to_cd)) == 0L) {
            MultiAssayExperiment::MultiAssayExperiment(
                experiments = to_exp,
                metadata = to_md
            )
        } else {
            MultiAssayExperiment::MultiAssayExperiment(
                experiments = to_exp,
                metadata = to_md,
                colData = to_cd
            )
        }
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
S7::method(as.list, MultiFactor) <- function(x, use.names = TRUE) {
    ifelse(
        use.names,
        yes = return(unfactor(x)),
        no = return(`names<-`(x@index, rownames(x)))
    )}
