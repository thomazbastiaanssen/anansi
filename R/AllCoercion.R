#' Coercion methods for AnansiWeb
#' @description
#' Coerce AnansiWeb to and from other object types.
#' @name coerceAnansi
#' @examples
#' #Create a random web
#' web <- randomWeb()
#'
#' #To list
#' as.list(web)
#'
#' #To MultiAssayExperiment
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
    out <- list(tableY = from@tableY, tableX = from@tableX,
                dictionary = from@dictionary, metadata = from@metadata)
    names(out)[c(1L, 2L)] <- names(from)
    out
})

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#'
setAs(from = "AnansiWeb", to = "MultiAssayExperiment", def = function(from) {

    tY  <- t(from@tableY)
    tX  <- t(from@tableX)
    to_exp <- ExperimentList(
        y = SummarizedExperiment(tY),
        x = SummarizedExperiment(tX)
    )
    names(to_exp) <- names(from)

    to_md  <- list(dictionary = from@dictionary)
    to_cd  <- from@metadata

    MultiAssayExperiment(
        experiments = to_exp,
        metadata = to_md,
        colData = to_cd
    )
})


