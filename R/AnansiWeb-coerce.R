#' @importFrom methods as
#' @export
#'
setAs(from = "anansiWeb", to = "list", def = function(from) {
    out <- list(from@tableY, from@tableX, dictionary = from@dictionary)
    names(out)[c(1L, 2L)] <- names(from)
    out
})

#' Coerce to MultiAssayExperiment
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
setAs(from = "anansiWeb", to = "MultiAssayExperiment", def = function(from) {
    weblist <- as(from, "list")
    experiments <- ExperimentList(
        y = SummarizedExperiment(t(weblist[[1L]])),
        x = SummarizedExperiment(t(weblist[[2L]]))
        )
    names(experiments) <- names(weblist)[c(1L, 2L)]
    metadata <- weblist[3L]

    MultiAssayExperiment(
        experiments = experiments,
        metadata = metadata)
    })
