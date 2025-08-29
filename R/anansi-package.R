#' anansi: Annotation-based Analysis of Specific Interactions.
#'
#' @description
#' The anansi package package provides tools to prepare and facilitate
#' integrative association analysis between the features of two data sets that
#' are known to interact.
#'
#' ## 1. Input for `anansi()` with [AnansiWeb()] and [MultiFactor()]
#' \itemize{
#'     \item [randomAnansi], [kegg_link()]: Generate example input
#'     \item [AnansiWeb], [MultiFactor]: Handle and manipulate
#'     input
#'     }
#'
#' ## 2. Output and cross-compatibility
#' \itemize{
#'     \item [`getGraph()`]: Compatibility with [igraph::igraph]
#'     \item [`plotAnansi()`]: Plot output in the style of [miaViz::miaViz]
#'     }
#'
#' ## 3. Vignettes
#' \itemize{
#'     \item \href{https://thomazbastiaanssen.github.io/anansi/articles/anansi.html}{1. Getting started with anansi}
#'     \item \href{https://thomazbastiaanssen.github.io/anansi/articles/adjacency_matrices.html}{2. Adjacency matrices}
#'     \item \href{https://thomazbastiaanssen.github.io/anansi/articles/association_testing.html}{3. Association testing}
#'     }
#'
#' @aliases anansi-package
#' @name anansi
#' @keywords internal
#' @docType package
"_PACKAGE"
