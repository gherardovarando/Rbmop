% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/1-base-functions.R
\name{summary.bmop}
\alias{summary.bmop}
\title{Summary of a bmop object}
\usage{
\method{summary}{bmop}(object, ...)
}
\arguments{
\item{object}{bmop object}

\item{...}{compatibility with \code{\link{summary}}}
}
\description{
Summary of a bmop object
}
\examples{
data(trees)
bmop<-bmop_fit(data=trees$Height)
summary(bmop)
}

