% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anansiCorTestByGroup.R
\name{anansiCorPvalue}
\alias{anansiCorPvalue}
\title{Compute r-statistics for each featureY-featureX pair in the dictionary.
Typically, the main \code{anansi()} function will run this for you.}
\usage{
anansiCorPvalue(web, group.bool, verbose)
}
\arguments{
\item{web}{An \code{AnansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.}

\item{group.bool}{A categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.}

\item{verbose}{A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.}
}
\value{
An \code{anansiTale} result object.
}
\description{
Compute r-statistics for each featureY-featureX pair in the dictionary.
Typically, the main \code{anansi()} function will run this for you.
}
\seealso{
\code{\link[=anansi]{anansi()}} \cr \code{\link[=anansiCorTestByGroup]{anansiCorTestByGroup()}}
}
