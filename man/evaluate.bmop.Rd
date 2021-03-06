% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/1-base-functions.R
\name{evaluate.bmop}
\alias{evaluate.bmop}
\title{Evaluation of a bmop object}
\usage{
evaluate.bmop(x, object, MIN = bmopPar()$MIN)
}
\arguments{
\item{x}{numeric, vector, matrix}

\item{object}{bmop object}

\item{MIN}{numeric}
}
\value{
Numeric value (values) of the computed bspline at point
 (points) \code{x}.
}
\description{
Evaluation of a bmop object
}
\details{
The deBoor algorithm is used, implemented in C.
}
\examples{
bmop<-bmop_fit(rnorm(100))
evaluate.bmop(0,bmop)
evaluate.bmop(c(-1,0,+1),bmop)
}
\references{
Carl de Boor,
On calculating with B-splines, \emph{Journal of Approximation Theory},
Volume 6, Issue 1, July 1972, Pages 50-62,
\url{http://www.sciencedirect.com/science/article/pii/0021904572900809}.
}

