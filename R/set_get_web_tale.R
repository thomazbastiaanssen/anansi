#' @noRd
#'
tell_F <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@f.values)
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_F))
  }
}

#' @noRd
#'
tell_T <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@t.values)
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_T))
  }
}

#' @noRd
#'
tell_P <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@p.values)
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_P))
  }
}


#' @noRd
#'
tell_e <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@estimates)
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_e))
  }
}


#' @noRd
#'
tell_df1 <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@df[1])
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_df1))
  }
}

#' @noRd
#'
tell_df2 <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@df[2])
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_df2))
  }
}

#' @noRd
#'
tell_dfr <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@df[3])
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_dfr))
  }
}

#' @noRd
#'
get_dict.double <- function(web) {
  `mode<-`(web@dictionary, "double")
}

#' @noRd
#' @importFrom Matrix as.matrix
get_dict.logical <- function(web) {
  `mode<-`(Matrix::as.matrix(web@dictionary), "logical")
}

#' @noRd
#' @importFrom Matrix as.matrix
get_dict <- function(web) {
  Matrix::as.matrix(web@dictionary)
}

#' @noRd
#'
get_tableX <- function(web) {
  web@tableX
}

#' @noRd
#'
get_tableY <- function(web) {
  web@tableY
}

