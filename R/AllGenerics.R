getEdgeList        <- S7::new_generic("getEdgeList", "x")
getGraph           <- S7::new_generic("getGraph", "x")
getFeaturePairs    <- S7::new_generic("getFeaturePairs", "x")

weaveWeb           <- S7::new_generic("getFeaturePairs", "x")
plotAnansi         <- S7::new_generic("plotAnansi", "x")


tableX      <- S7::new_generic("tableX", "x")
tableY      <- S7::new_generic("tableY", "x")
`tableX<-`  <- S7::new_generic("tableX<-", "x")
`tableY<-`  <- S7::new_generic("tableY<-", "x")
dictionary     <- S7::new_generic("dictionary", "x")
`dictionary<-` <- S7::new_generic("dictionary<-", "x")

metadata       <- S7::new_external_generic("S4Vectors", "metadata", "x")
which          <- S7::new_external_generic("BiocGenerics", "which", "x")
mapply         <- S7::new_external_generic("BiocGenerics", "mapply", "...")

as.data.frame  <- S7::new_external_generic("base", "as.data.frame", "x")


# Externals
unfactor     <- S7::new_external_generic("S4Vectors", "unfactor", "x")
# droplevels   <- S7::new_external_generic("base","droplevels", "x")
# levels       <- S7::new_external_generic("base", "levels", "x")
# as.list      <- S7::new_external_generic("base", "as.list", "x")
# show         <- S7::new_external_generic("methods", "show", "object")
# names        <- S7::new_external_generic("base", "names", "x")
# dim          <- S7::new_external_generic("base", "dim", "x")
# dimnames     <- S7::new_external_generic("base", "dimnames", "x")





