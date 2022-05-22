#' Make a progress bar. Lifted from the excellent ALDEx2 and exprso libraries.
#' @description plot a progress bar.
#' @param i current state
#' @param k maximum
#' @param numTicks How many bars to plot
#' @return a progress bar
#'
progress <- function(i, k, numTicks){

  if(i == 1) numTicks <- 0

  if(numTicks == 0) cat("|-")

  while(i > numTicks*(k/40)){

    cat("-")
    if(numTicks == 10) cat("(25%)")
    if(numTicks == 20) cat("(50%)")
    if(numTicks == 30) cat("(75%)")
    numTicks <- numTicks + 1
  }

  if(i == k) cat("-|\n")

  return(numTicks)
}
