#' @title Coercion functions for anansi
#' @name as.list.AnansiWeb
#' @rdname anansi-coercion
#' @method as.list AnansiWeb
#' @returns An object of the desired class.
#' @param x input object
#' @param use.names `Logical scalar`, whether output list should contain
#'     character (Default) or integer data frame. If `FALSE`, returns
#'     `unfactor(x)`.
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
#'
#' # MultiFactor
#' x <- randomMultiFactor()
#' as.list(x, use.names = TRUE)
#'
S7::method(convert, list(AnansiWeb, S7::class_list)) <-
               function(from, to) as.list.AnansiWeb(x = from)

#' @importFrom S7 convert
#' @method as.data.frame AnansiWeb
#' @rdname anansi-coercion
#' @name as.data.frame.AnansiWeb
#'
S7::method(convert, list(AnansiWeb, S7::class_data.frame)) <-
    function(from, to) as.data.frame.AnansiWeb(x = from)

#' @name as.list.MultiFactor
#' @rdname anansi-coercion
#' @importFrom S7 convert
#'
S7::method(convert, list(MultiFactor, S7::class_list)) <-
    function(from, to) as.list(from)

#' @export
S7::method(as.list, MultiFactor) <- function(x, ..., use.names = TRUE) {
    ifelse(
        use.names,
        yes = return(S4Vectors::unfactor(x)),
        no = return(`names<-`(x@index, rownames(x)))
    )}

#' @export
S7::method(as.list, AnansiWeb) <- function(x, ...) {
    out <- S7::props(x)
    names(out)[c(1L, 2L)] <- names(x)
    out
}

#' @export
S7::method(as.data.frame, AnansiWeb) <- function(x, row.names, optional, ...) {
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
#' @usage NULL
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


