\name{Rbmop-package}
\alias{Rbmop-package}
\alias{Rbmop}
\docType{package}
\title{
Handling and estimation of B-spline MoPs densities and conditional densities
(bmop)
}
\description{
Handling and estimation of B-spline MoPs densities and conditional densities 
(bmop)
}
\details{
\tabular{ll}{
Package: \tab Rbmop\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2015-06-10\cr
License: \tab GPL-2 \cr
}
}
\author{
Gherardo Varando
Maintainer: Gherardo Varando <gherardo.varando@gmail.com>
}
\references{
...
}
\keyword{ package }
\seealso{
}
\examples{
data(trees)
bmop<-bmop_fit(data=trees$Height)
summary(bmop)
envelope_plot(n=50,N=50,rtrue=rexp,dtrue=dexp)
require("MASS")
X<-rnorm(100)
Y<-rnorm(100,mean=X)
data<-data.frame(Y,X)
condbmop<-bmop_fit(data=data,conditional=TRUE)
plot(condbmop,N=20,persp=TRUE)
}
