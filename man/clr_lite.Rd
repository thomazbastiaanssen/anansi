% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clr_lite.R
\name{clr_lite}
\alias{clr_lite}
\alias{clr_c}
\alias{clr_unif}
\alias{clr_logunif}
\title{Impute zeroes and perform a centered log-ratio (CLR) transformation}
\usage{
clr_lite(counts, samples_are = "cols", method = "logunif", replicates = 1000)

clr_c(counts, samples_are = "cols")

clr_unif(counts, samples_are = "cols", replicates = 1000)

clr_logunif(counts, samples_are = "cols", replicates = 1000)
}
\arguments{
\item{counts}{A compositional count table.}

\item{samples_are}{Either "cols" or "rows". Default is "cols". Denotes whether the columns or rows depict individual samples.}

\item{method}{The method for zero imputation. One of \code{"logunif"}, \code{"unif"} or \code{"const"}.
\code{'logunif'} samples small numbers from a log-uniform distribution, whereas \code{'unif'} samples from a uniform one. On the other hand, \code{"const"} simply replaces zeroes with \code{0.65 * [the lowest value]}.}

\item{replicates}{An integer. For the two random sampling methods, if this is larger than 1, every zero will be imputed that many times. The median of the CLR of all those replicates will be returned. If \code{method} is set to \code{"const"}, replicates will be automatically set to 1 as no random numbers are generated.}
}
\value{
A CLR-transformed count table.
}
\description{
Microbiome data is compositional. When compositional data is examined using non-compositional methods, many problems arise.
Performing a centered log-ratio transformation is a reasonable way to address these problems reasonably well.
\cr \cr
A major problem with this approach is that microbiome data typically contains lots of zeroes and the logarithm of zero is undefined.
Here, we implemented a few methods discussed by Lubbe \emph{et al.} 2021 to replace zeroes with non-zero values in such a way that the structure of the data remains reasonably well preserved.
\cr \cr
Some of these methods (namely 'logunif' and 'runif') involve imputing small values between 0 and the lowest non-zero value in the dataset.
For these methods, we have implemented a resampling approach in order to stabilize the inter-run variability. See \code{method} for more information.
}
\section{Functions}{

\code{clr_c:}
 A wrapper for \code{clr_lite(counts, method = "const", replicates = 1)}.


\code{clr_unif:}
 A wrapper for \code{clr_lite(counts, method = "unif")}.


\code{clr_logunif:}
 A wrapper for \code{clr_lite(counts, method = "logunif")}.
}

\references{
Sugnet Lubbe, Peter Filzmoser, Matthias Templ (2021)
\emph{Comparison of zero replacement strategies for compositional data with large numbers of zeros.}
doi:https://doi.org/10.1016/j.chemolab.2021.104248
}
