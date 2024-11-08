#' @noRd
#'
get_dict.double <- function(web)
  `mode<-`(web@dictionary, "double")

#' @noRd
#'
get_dict.logical <- function(web)
  `mode<-`(web@dictionary, "logical")

#' @noRd
#'
get_dict <- function(web)
  web@dictionary

#' @noRd
#'
get_tableX <- function(web)
  web@tableX

#' @noRd
#'
get_tableY <- function(web)
  web@tableY
