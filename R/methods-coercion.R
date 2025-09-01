#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases as.list,anansi::AnansiWeb-method
#' @usage NULL
#' @method as.list AnansiWeb
#'
S7::method(as.list, AnansiWeb)  <- function(x, ...) {
    out <- S7::props(x)
    names(out)[c(1L, 2L)] <- names(x)
    out
}

#' @name AnansiWeb
#' @rdname AnansiWeb
#' @method as.data.frame `anansi::AnansiWeb`
#' @usage NULL
#'
`as.data.frame.anansi::AnansiWeb` <- function(x, ...) {
    cbind(
        x@tableY,
        x@tableX,
        x@metadata
    )
}

#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases as.data.frame.anansi::AnansiWeb
#' @importFrom S4Vectors as.data.frame
#' @importFrom S7 convert
#' @export
#'
S7::method(convert, list(AnansiWeb, S7::class_data.frame)) <-
    function(from, to) `as.data.frame.anansi::AnansiWeb`(from)

#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases as.list.anansi::AnansiWeb
#' @importFrom S4Vectors as.data.frame
#' @importFrom S7 convert
#' @export
#'
S7::method(convert, list(MultiFactor, S7::class_list)) <-
               function(from, to) `as.list.anansi::AnansiWeb`(x = from)

#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases as.MAE asMAE as.MultiAssayExperiment asMultiAssayExperiment
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
#' @method as.list `anansi::MultiFactor`
#'
`as.list.anansi::MultiFactor` <- function(x, ..., use.names = TRUE) {
    ifelse(
        use.names,
        yes = return(S4Vectors::unfactor(x)),
        no = return(`names<-`(x@index, rownames(x)))
    )}

