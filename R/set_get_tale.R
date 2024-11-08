#' @noRd
#'
tell_F <- function(tale)
  tale@F.values

#' @noRd
#'
tell_P <- function(tale)
  tale@p.values

#' @noRd
#'
tell_e <- function(tale)
  tale@estimates

#' @noRd
#'
tell_df1 <- function(tale)
  tale@df[1]

#' @noRd
#'
tell_df2 <- function(tale)
  tale@df[2]


#' @noRd
#'
tell_dfr <- function(tale)
  tale@df[3]
