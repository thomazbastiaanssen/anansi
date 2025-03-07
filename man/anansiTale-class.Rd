% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{anansiTale-class}
\alias{anansiTale-class}
\title{An S4 class to contain all \code{anansi} stats results so that they can
easily be extracted.}
\description{
\code{anansiTale} is the main container that will hold your
stats output data coming out of the \code{anansi} pipeline.
}
\section{Slots}{

\describe{
\item{\code{subject}}{A character that describes the data that was queried.}

\item{\code{type}}{A character that describes type of parameter contained in the
\code{estimates} slot. For example r.values for correlations or r.squared
for models.}

\item{\code{df}}{a vector of length 2, containing df1 and df2 corresponding to the
F-ratio considered.}

\item{\code{estimates}}{A matrix containing the estimates for the parameters named in
the \code{type} slot.}

\item{\code{f.values}}{A matrix containing the f-values, for least-squares.}

\item{\code{t.values}}{A matrix containing the t-values, for correlations.}

\item{\code{p.values}}{A matrix containing the p.values for the parameters named in
the \code{type} slot.}
}}

