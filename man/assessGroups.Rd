% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anansi.R
\name{assessGroups}
\alias{assessGroups}
\title{Investigate validity of the groups argument}
\usage{
assessGroups(web, groups, diff_cor = diff_cor)
}
\arguments{
\item{web}{web An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.}

\item{groups}{A categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.}

\item{diff_cor}{A boolean. Toggles whether to compute differential correlations. Default is \code{TRUE}.}
}
\value{
a list including a modified \code{groups} and \code{diff_cor} argument.
}
\description{
checks whether \code{groups} is missing, the correct length and suitable for correlations per groups and differential correlations.
If something looks off, \code{assessGroups} will do its best to salvage it and let you know something's up.
This is a helper function called by \code{anansi}.
}
