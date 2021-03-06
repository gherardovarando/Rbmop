% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/2-estimation-functions.R
\name{search_bmop}
\alias{search_bmop}
\title{Greedy penalized log-likelihood search}
\usage{
search_bmop(data, conditional = F, k = Rbmop::bmopPar()$k,
  corrected = FALSE, knotsMethod = Rbmop::bmopPar()$knotsMethod, ...)
}
\arguments{
\item{data}{data.frame, matrix or vector, the variables must be in
the right order (the columns of data)}

\item{conditional}{logical}

\item{k}{positive number or \code{"BIC"}}

\item{corrected}{logical}

\item{knotsMethod}{the method to use in knots generation}

\item{...}{additional parameters}
}
\value{
A bmop object, the aproximations of \eqn{f}.
}
\description{
Aproximation of a density \eqn{f(x_1,\ldots,x_n)} or conditional density
}
\examples{
data<-rnorm(100)
bmopS<-search_bmop(data=data)
plot(bmopS)
}

