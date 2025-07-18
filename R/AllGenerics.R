#' @export
getEdgeList        <- S7::new_generic("getEdgeList", "x")

#' @export
getGraph           <- S7::new_generic("getGraph", "x")

#' @export
getFeaturePairs    <- S7::new_generic("getFeaturePairs", "x")

#' @export
weaveWeb           <- S7::new_generic("weaveWeb",   "x")

#' @export
plotAnansi         <- S7::new_generic("plotAnansi", "x")

#' @export
tableX         <- S7::new_generic("tableX",         "x")

#' @export
tableY         <- S7::new_generic("tableY",         "x")

#' @export
`tableX<-`     <- S7::new_generic("tableX<-",       "x")

#' @export
`tableY<-`     <- S7::new_generic("tableY<-",       "x")

#' @export
dictionary     <- S7::new_generic("dictionary",     "x")

#' @export
`dictionary<-` <- S7::new_generic("dictionary<-",   "x")

#' @export
metadata       <- S7::new_external_generic("S4Vectors", "metadata", "x")

#' @export
`metadata<-`   <- S7::new_external_generic("S4Vectors", "metadata<-", "x")


# Workarounds
#' @export
which <- S7::new_generic("which", "x")

# For every other class (base, S3, S4, S7), run BiocGenerics
#' @export
S7::method(which, S7::class_any) <- function(x) {
    BiocGenerics::which(x)
}


#' Workarounds
#' @export
pairwiseApply <- S7::new_generic("pairwiseApply", "X")

#as.data.frame  <- S7::new_external_generic("BiocGenerics", "as.data.frame", "x")
#' @export
as.data.frame  <- S7::new_generic("as.data.frame", "x")

#' @export
S7::method(as.data.frame, S7::class_any) <- function(x, ...) {
    BiocGenerics::as.data.frame(x, ...)
}

# Externals
# unfactor     <- S7::new_external_generic("S4Vectors", "unfactor", "x")
#' @export
unfactor     <- S7::new_generic("unfactor", "x")

# droplevels   <- S7::new_external_generic("base","droplevels", "x")
# levels       <- S7::new_external_generic("base", "levels", "x")
# as.list      <- S7::new_external_generic("base", "as.list", "x")
# show         <- S7::new_external_generic("methods", "show", "object")
# names        <- S7::new_external_generic("base", "names", "x")
# dim          <- S7::new_external_generic("base", "dim", "x")
# dimnames     <- S7::new_external_generic("base", "dimnames", "x")
