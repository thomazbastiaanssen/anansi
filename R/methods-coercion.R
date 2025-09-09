#' @title Coercion functions for anansi
#' @name anansi-coercion
#' @rdname anansi-coercion
#' @returns An object of the desired class.
#' @param x input object
#' @param use.names `Logical scalar`, whether output list should contain
#'     character (Default) or integer data frame. If `FALSE`, returns
#'     `unfactor(x)`.
#' @param row.names,optional Not used. See ?base::data.frame
#' @param ... additional arguments (currently not used).
#' @seealso [unfactor()]
#' @examples
#' # AnansiWeb
#' x <- randomWeb(36)
#'
#' as.list(x)
#' as.data.frame(x)
#'
#' # AnansiWeb to MultiAssayExperiment
#' asMAE(x)
#' asTSE(x)
#'
#' # MultiFactor
#' x <- randomMultiFactor()
#' as.list(x, use.names = TRUE)
#'
NULL

#' @export
S7::method(convert, list(AnansiWeb, S7::class_list)) <-
    function(from, to) `as.list.anansi::AnansiWeb`(x = from)

#' @export
#' @rdname anansi-coercion
`as.list.anansi::AnansiWeb` <- function(x, ...) {
    out <- S7::props(x)
    names(out)[c(1L, 2L)] <- names(x)
    out
}

#' @importFrom S7 convert
#' @export
#'
S7::method(convert, list(AnansiWeb, S7::class_data.frame)) <-
    function(from, to) as.data.frame.AnansiWeb(x = from)

#' @importFrom S7 convert
#' @export
#'
S7::method(convert, list(MultiFactor, S7::class_list)) <-
    function(from, to) `as.list.anansi::MultiFactor`(x = from)

#' @export
S7::method(as.list, MultiFactor) <-
    function(x, ..., use.names = TRUE) {
        `as.list.anansi::MultiFactor`(
            x, ..., use.names
        )
    }

#' @export
#' @importFrom S4Vectors unfactor
#' @rdname anansi-coercion
`as.list.anansi::MultiFactor` <- function(x, ..., use.names = TRUE) {
    ifelse(
        use.names,
        yes = return(S4Vectors::unfactor(x)),
        no = return(`names<-`(x@index, rownames(x)))
    )
}

#' @export
S7::method(as.list, AnansiWeb) <- function(x, ...) {
    `as.list.anansi::AnansiWeb`(x)
}

#' @export
S7::method(as.data.frame, AnansiWeb) <- function(x, row.names, optional, ...) {
    `as.data.frame.anansi::AnansiWeb`(x)
}

#' @export
#' @rdname anansi-coercion
`as.data.frame.anansi::AnansiWeb` <- function(x, row.names, optional, ...) {
    cbind(
        x@tableY,
        x@tableX,
        x@metadata
    )
}

#' @name asMAE
#' @rdname anansi-coercion
#' @aliases as.MAE as.MultiAssayExperiment asMultiAssayExperiment
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#'
asMAE <- function(x) {
    tY <- t(x@tableY)
    tX <- t(x@tableX)
    to_exp <- ExperimentList(
        y = SummarizedExperiment::SummarizedExperiment(tY),
        x = SummarizedExperiment::SummarizedExperiment(tX)
    )
    names(to_exp) <- names(x)

    to_md <- list(dictionary = x@dictionary)
    to_cd <- x@metadata

    if (prod(dim.data.frame(to_cd)) == 0L) {
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

#' @name asTSE
#' @rdname anansi-coercion
#' @aliases as.TSE as.TreeSummarizedExperiment asTreeSummarizedExperiment
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @export
#'
asTSE <- function(x) {
    tY <- t(x@tableY)
    tX <- t(x@tableX)
    to_exp <- list(
        y = tY
    )
    names(to_exp) <- names(x)[1]
    to_alt <- list(
        x = TreeSummarizedExperiment::TreeSummarizedExperiment(tX)
    )
    names(to_alt) <- names(x)[2]

    to_md <- list(dictionary = x@dictionary)
    to_cd <- x@metadata

    if (prod(dim.data.frame(to_cd)) == 0L) {
        TreeSummarizedExperiment::TreeSummarizedExperiment(
            assays = to_exp,
            altExps = to_alt,
            metadata = to_md
        )
    } else {
        TreeSummarizedExperiment::TreeSummarizedExperiment(
            assays = to_exp,
            altExps = to_alt,
            metadata = to_md,
            colData = to_cd
        )
    }
}
