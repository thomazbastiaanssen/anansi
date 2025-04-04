#' @rdname getWeb
#' @export
#'
#' @param x a [MultiAssayExperiment::MultiAssayExperiment()].
#'
#' @param tableY,tableX
#' `Character scalar` or `numeric scalar`. Selects experiment
#' corresponding to `tableY` and `tableX` from `experiments(x)`
#' of `MultiAssayExperiment` object by name or index, name is recommended.
#' (Default slots: `Y = 1`, `X = 2`).
#'
#' @param link One of the following:
#' \itemize{
#'  \item Character scalar with value "none"
#'  \item data.frame with two columns
#'  \item list with two such data.frames
#' }
#'
#' @param force_new `boolean` If x already has a dictionary `Matrix`
#' in metadata, ignore it and generate a new object anyway? (Default: FALSE).
#' @param ... additional parameters passed to [AnansiWeb()].
#' @param tableY,tableX
#' `Character scalar` or `numeric scalar`. Selects experiment
#' corresponding to `tableY` and `tableX` from `experiments(x)`
#' of `MultiAssayExperiment` object by name or index, name is recommended.
#' (Default slots: `Y = 1L`, `X = 2L`).
#' @param typeY,typeX `Character scalar` or `numeric scalar`. Selects assay from
#' experiments to `tableY` and `tableX` from `experiments(x)`.
#' (Default: `1L` - the first assay in that experiment).
#' @param experiment1,experiment2 synonymous args to `tableY,tableX` for
#'     compatibility with `mia` argument style.
#' @param assay.type1,assay.type2 synonymous args to `typeY,typeX` for
#'     compatibility with `mia` argument style.
#' @returns an `AnansiWeb` object, with sparse binary biadjacency matrix
#' with features from `y` as rows and features from `x` as columns in
#' `dictionary` slot. If x already contains a dictionary in metadata, use
#' that one, unless `force_new = TRUE`.
#'
#' @importFrom MultiAssayExperiment MultiAssayExperiment metadata experiments
#' @importFrom SummarizedExperiment assay colData
#' @importClassesFrom Matrix Matrix
#'
#' @details
#' This wrapper of [weaveWeb()] allows to generate an
#' `AnansiWeb` S4 object directly from objects of class
#' [MultiAssayExperiment::MultiAssayExperiment()]
#' . First, the assays specified by `assay.typeY` and `assay.typeX`
#' are passed to [AnansiWeb()] to build an AnansiWeb object.
#'
#' @examples
#' # Make a random anansiWeb
#' web <- randomWeb()
#'
#' # Combine experiments into MultiAssayExperiment object
#' mae <- as(web, "MultiAssayExperiment")
#'
#' # Back to AnansiWeb
#' outWeb <- getWeb(mae, tableY = "y", tableX = "x")
#'
setMethod("getWeb",
          signature = c(x = "MultiAssayExperiment"),
          function(
        x, link = NULL, ...,
        tableY = NULL, tableX = NULL,
        typeY = NULL, typeX = NULL,
        force_new = FALSE,
        experiment1 = NULL, experiment2 = NULL,
        assay.type1 = NULL, assay.type2 = NULL
    ) {
      y_ids <- .test_coherent(tableY, experiment1, typeY, assay.type1)
      x_ids <- .test_coherent(tableX, experiment2, typeX, assay.type2)
      tableY <- y_ids[[1L]]
      tableX <- x_ids[[1L]]
    # Check experiments
    mia:::.test_experiment_of_mae(x, tableY)
    mia:::.test_experiment_of_mae(x, tableX)
    y_exp <- names(experiments(x)[tableY])
    x_exp <- names(experiments(x)[tableX])

    # Extract assays
    tY <- t(assay(experiments(x)[tableY], y_ids[[2L]]))
    tX <- t(assay(experiments(x)[tableX], x_ids[[2L]]))

    # Check if x already contains a dictionary
    if (!force_new) {
        m <- metadata(x)

        if (is.null(link)) {
            d <- "dictionary"
        }
        if (valid_selection(link, m)) {
            d <- link
        }
        if (d %in% names(m)) {
            return(AnansiWeb(
                tableX = tX, tableY = tY, dictionary = m[[d]],
                metadata = list(metadata = as.data.frame(colData(x))),
                ...
            ))
        }
    }
    # Generate web object
    weaveWeb.default(
        x = x_exp, y = y_exp, link = link,
        tableX = tX, tableY = tY,
        metadata = list(metadata = as.data.frame(colData(x))),
        ...
    )
    }
)

#' TRUE if i can select in x
#' @noRd
#' @param i `Character or numeric scalar`. Index to check.
#' @param x object to check i in
#' @returns TRUE if i selects in x, FALSE otherwise
#'
valid_selection <- function(i, x) {
    # Need to be length 1.
    if (length(i) != 1L) FALSE

    # If numeric, needs to be within element length of x
    if (is.numeric(i)) i <= length(x)
    # If character, x needs to be named and i needs to be within those names
    if (is.character(i)) {
        if (is.null(names(x))) FALSE

        !is.na(match(i, names(x)))
    }
    # If that didn't work, invalid selection. return FALSE.
    FALSE
}

#' @noRd
.test_coherent <- function(tab, exp, tab.type, ass.type) {
    e_out <- unique(c(tab, exp))
    if(length(e_out) == 0L) {e_out <- 1L}
    stopifnot(
        "args 'tableY,X' cannot contradict args 'experiment1,2'." =
                  length(e_out) == 1L
        )

    t_out <- unique(c(tab.type, ass.type))
    if(length(t_out) == 0L) {t_out <- 1L}
    stopifnot(
        "args 'typeY,X.' cannot contradict args 'assay.type1,2'." =
                  length(t_out) == 1L
        )

    return(list(e_out, t_out))
}
