% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ParseOutput.R
\name{spinToLong}
\alias{spinToLong}
\title{Take anansi output and wrangle it all to a long data.frame format.}
\usage{
spinToLong(
  anansi_output,
  prune = T,
  translate = F,
  Y_translation = NULL,
  X_translation = NULL
)
}
\arguments{
\item{anansi_output}{an \code{anansiYarn} object. The output of the main anansi function.}

\item{prune}{Boolean, default is TRUE. Toggles whether to take out the non-canonical associations.}

\item{translate}{Boolean, default is FALSE Toggles whether to translate the names of the features in tableX and tableY to human readable names.}

\item{Y_translation}{data.frame, a lookup table with featureY names as the first column and human readable names as the second. See \code{cpd_translation} for an example file.}

\item{X_translation}{data.frame, a lookup table with featureX names as the first column and human readable names as the second. See \code{KO_translation} for an example file.}
}
\value{
a long format data.frame intended to be compatible with \code{ggplot2}
}
\description{
Take anansi output and wrangle it all to a long data.frame format.
}
