#' AnansiWeb S4 container class
#' @name AnansiWeb
#' @rdname AnansiWeb
#' @aliases AnansiWeb-class
#' @slot tableY,tableX Two `matrix` objects of measurements, data. Rows are
#'     samples and columns are features. Access with `tableY()` and `tableX()`.
#' @slot dictionary `Matrix`, binary adjacency matrix. Optionally sparse.
#'     Typically generated using the `weaveWeb()` function. Access with
#'     `dictionary()`.
#' @slot metadata Optional `data.frame` of sample metadata. Access with
#'     `metadata()`.
#' @importClassesFrom Matrix Matrix
#' @importClassesFrom S4Vectors Annotated
#' @importMethodsFrom S4Vectors metadata
#'
setClass(
    "AnansiWeb",
    contains = "Annotated",
    slots = c(
        tableY = "matrix",
        tableX = "matrix",
        dictionary = "Matrix",
        metadata = "list"
    )
)

#' is valid AnansiWeb?
#' @noRd
#' @description
#' returns TRUE if input is in the right format to be an AnansiWeb object
#' @param object
#' `any` object, but not much will happen unless the object's class has a
#' formal definition.
#' @importFrom methods validObject
#' @returns `TRUE` if passes, character vector otherwise.
#'
setValidity("AnansiWeb", method = function(object) {
    ifelse(
        test = validWeb(object),
        yes = TRUE,
        no = "object is not in a valid format."
    )
})

#' MultiFactor S4 container class
#' @rdname MultiFactor
#' @name MultiFactor
#' @aliases MultiFactor-class
#' @slot index Named `list` of named integer data frames of at least two columns
#'     each. The column names correspond to names in the `levels` slot. Similar
#'     to `factor`s, the integers in those columns correspond to the characters
#'     in that level. Accessed through regular list methods (e.g., `[`, `[[`).
#' @slot levels `Named list of character vectors`. Accessed through `levels(x)`
#' @slot map `(sparse)Matrix` specifying which elements contain which levels.
#'     Accesses through `dictionary(x)`.
#' @importClassesFrom Matrix Matrix
#' @export
#'
setClass(
    "MultiFactor",
    slots = c(
        index = "list",
        levels = "list",
        map = "Matrix"
    )
)

#' is valid MultiFactor?
#' @noRd
#' @description
#' returns TRUE if input is in the right format to be an MultiFactor object
#' @param object
#' `any` object, but not much will happen unless the object's class has a
#' formal definition.
#' @importFrom methods validObject
#' @returns `TRUE` if passes, character vector otherwise.
#'
setValidity("MultiFactor", method = function(object) {
    ifelse(
        test = validMultiFactor(object),
        yes = TRUE,
        no = "object is not in a valid format."
    )
})


#' Is this a data.frame with exactly two columns that are named?
#' @noRd
validMultiFactor <- function(x) {
    levels_valid <- validLevels(levels(x))
    x <- x@index
    values_valid <- vapply(x, validIntLinkDF, NA, USE.NAMES = FALSE)
    no_missing <- !any(vapply(x, anyNA, NA, USE.NAMES = FALSE))

    if (!isTRUE(levels_valid)) {
        message("Levels are not structured correctly. ")
    }
    if (!isTRUE(all(values_valid))) {
        message(
            "List content in positions ",
            paste(which(!isTRUE(values_valid)), collapse = ", "),
            " not structured correctly. "
        )
    }
    if (!no_missing) {
        message("Missing values are not allowed.")
    }

    if (
        all(
            levels_valid,
            isTRUE(values_valid),
            isTRUE(all(no_missing))
        )
    ) {
        return(TRUE)
    }
}

#' An S4 class to contain all `anansi` stats results so that they can
#' easily be extracted.
#'
#' @slot subject A character that describes the data that was queried.
#' @slot type A character that describes type of parameter contained in the
#'     `estimates` slot. For example r.values for correlations or r.squared
#'     for models.
#' @slot df a vector of length 2, containing df1 and df2 corresponding to the
#'     F-ratio considered.
#' @slot estimates A matrix containing the estimates for the parameters named in
#'     the `type` slot.
#' @slot f.values A matrix containing the f-values, for least-squares.
#' @slot t.values A matrix containing the t-values, for correlations.
#' @slot p.values A matrix containing the p.values for the parameters named in
#'     the `type` slot.
#' @description `anansiTale` is the main container that will hold your
#'     stats output data coming out of the `anansi` pipeline.
#'
setClass(
    "anansiTale",
    slots = c(
        subject = "character",
        type = "character",
        df = "numeric",
        estimates = "matrix",
        f.values = "matrix",
        t.values = "matrix",
        p.values = "matrix"
    )
)
