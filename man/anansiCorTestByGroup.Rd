% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anansiCorTestByGroup.R
\name{anansiCorTestByGroup}
\alias{anansiCorTestByGroup}
\title{Run correlations for all interacting metabolites and functions.}
\usage{
anansiCorTestByGroup(web, group.vec, verbose = TRUE)
}
\arguments{
\item{web}{An \code{AnansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.}

\item{group.vec}{A character vector denoting group membership. Typically a state or treatment score.}

\item{verbose}{A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.}
}
\value{
a list of \code{anansiTale} result objects, one for the total dataset and per group if applicable.
}
\description{
If the \code{groups} argument is suitable, will also run correlation analysis per group. Typically, the main \code{anansi()} function will run this for you.
}
\seealso{
\code{\link[=anansi]{anansi()}}
}
