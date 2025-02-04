#' @noRd
#'
tell_F <- function(tale) {
  if (is(tale, "anansiTale")) {
    return(tale@F.values)
  }
  if (is.list(tale)) {
    return(lapply(tale, tell_F))
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
#'
get_dict.logical <- function(web) {
  `mode<-`(web@dictionary, "logical")
}

#' @noRd
#'
get_dict <- function(web) {
  web@dictionary
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

#' @noRd
#'
yarn.e <- function(yarn) {
  yarn@input@error.term
}

#' @noRd
#'
yarn.int <- function(yarn) {
  yarn@input@int.terms
}

#' @noRd
#'
yarn.grp <- function(yarn) {
  yarn@input@groups
}

#' @noRd
#'
yarn.f <- function(yarn) {
  yarn@input@lm.formula
}

#' @noRd
#'
yarn.web <- function(yarn) {
  yarn@input@web
}

#' @noRd
#'
yarn.dic <- function(yarn) {
  get_dict(yarn.web(yarn))
}

#' @noRd
#'
yarn.dic.double <- function(yarn) {
  get_dict.double(yarn.web(yarn))
}

#' @noRd
#'
yarn.dic.logical <- function(yarn) {
  get_dict.logical(yarn.web(yarn))
}


#' @noRd
#'
yarn.tX <- function(yarn) {
  get_tableX(yarn.web(yarn))
}

#' @noRd
#'
yarn.tY <- function(yarn) {
  get_tableY(yarn.web(yarn))
}




