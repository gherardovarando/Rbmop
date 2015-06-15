#' Rbmop: A package for handling and estimating densities and conditional
#' densities with B-spline.
#'
#' The main functions of Rbmop package are the following:
#' \itemize{
#' \item \code{\link{bmop_fit}} for estimation 
#' of density or conditional density.
#' \item \code{\link{bmopPar}} function for changing the parameter of the 
#' package.
#' \item \code{\link{compare.bmop}} and \code{\link{comaprison_plot}} comparison
#' functions.
#' \item Moreover the package provide base functions as \code{plot, print} and
#' conversion methods \code{\link{as.bmop}}.
#' }
#' 
#' @author Gherardo Varando \email{gherardo.varando@@gmail.com}, 
#' Concha Bielza and Pedro Larrañaga
#' @examples
#'data(trees)
#'bmop<-bmop_fit(data=trees$Height)
#'summary(bmop)
#'##################
#'
#'envelope_plot(n=50,N=50,rtrue=rexp,dtrue=dexp)
#'##################
#'## 
#'require("MASS")
#'X<-rnorm(100)
#'Y<-rnorm(100,mean=X)
#'data<-data.frame(Y,X)
#'condbmop<-bmop_fit(data=data,conditional=TRUE)
#'plot(condbmop,N=20,persp=TRUE)
#'
#'###################
#'## Trimodal-not differentiable density estimation via histogram
#'f<-function(x){return( (dnorm(x,mean=6)+dnorm(x,mean=-3)+dgamma(x,shape=2))/3 )}
#'D<-c(rnorm(10000,mean=-3),rgamma(10000,shape=2),rnorm(10000,mean=+6))
#'bmopPar(mle=T)
#'bmop<-bmop_fit(hist(D,nclass.FD))
#'comparison_plot(bmop,f,ylim=NULL)
#'
#' @docType package
#' @name Rbmop
NULL
#> NULL