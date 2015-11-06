#' Rbmop: A package for handling and estimating densities and conditional
#' densities with B-spline.
#'
#' The main functions of Rbmop package are the following:
#' \itemize{
#' \item \code{\link{bmop_fit}} for estimation 
#' of density or conditional density.
#' \item \code{\link{bmopPar}} function for changing the parameter of the 
#' package.
#' \item \code{\link{squareError.bmop}} and \code{\link{comparison_plot}} comparison
#' functions.
#' \item Moreover the package provide base functions as \code{plot, print} and
#' conversion methods \code{\link{as.function.bmop}}.
#' }
#' 
#' @author Gherardo Varando \email{gherardo.varando@@gmail.com}
#' @examples
#'data(trees)
#'bmop<-bmop_fit(data=trees$Height)
#'summary(bmop)
#'
#'##################
#'## 
#'X<-rnorm(100)
#'Y<-rnorm(100,mean=X)
#'data<-data.frame(Y,X)
#'condbmop<-bmop_fit(data=data,conditional=TRUE)
#'plot(condbmop,N=20,persp=TRUE)
#'
#'###################
#'## Trimodal-not differentiable density estimation via histogram
#'f<-function(x){return( (dnorm(x,mean=6)+
#'                        dnorm(x,mean=-3)+
#'                        dgamma(x,shape=2))/3 )}
#'D<-c(rnorm(1000,mean=-3),rgamma(1000,shape=2),rnorm(1000,mean=+6))
#'bmopPar(mle=TRUE)
#'bmop<-bmop_fit(hist(D,nclass.FD))
#'comparison_plot(bmop,f,ylim=NULL)
#'bmopPar(mle=FALSE)
#'
#'########################
#'## Using the examples functions
#'
#'Ex<-ex_bmop_gaussian2Mixture()
#'points(density(Ex$data),type="l",col="blue",lwd=3)
#'
#'####
#'
#' Ex<-ex_bmop_gaussianBetaGamma(N = 10000,m1 = 8)
#' points(density(Ex$data),type="l",col="blue",lwd=3)
#' squareError.bmop(Ex$bmop,Ex$true)
#' 
#' #######################
#' ## Envelopes
#' 
#' #Gaussian 
#' results<-envelope_plot(n=100,N=100,dtrue=dnorm,rtrue=rnorm)
#' print(results)
#' # Maximum Likelihood
#' bmopPar(mle=TRUE)
#' resultsMLE<-envelope_plot(n=100,N=100,dtrue=dnorm,rtrue=rnorm)
#' print(resultsMLE)
#' bmopPar(mle=FALSE)
#' 
#' #Exponential
#' 
#' results<-envelope_plot(n=100,N=1000,dtrue=dexp,rtrue=rexp)
#' print(results)
#' # Maximum Likelihood
#' bmopPar(mle=TRUE)
#' resultsMLE<-envelope_plot(n=100,N=100,dtrue=dexp,rtrue=rexp)
#' print(resultsMLE)
#' bmopPar(mle=FALSE)
#' 
#' #Gamma
#' 
#' envelope_plot(1000,1000,rtrue = rgamma,dtrue = dgamma,shape=3)
#' 
#' ############################################
#' ## Study of Errors and KL
#' 
#' Ns<-seq(from = 100,to=100000,by = 1000)
#' dim(Ns)<-c(length(Ns),1)
#' res<-as.data.frame(t(apply(Ns,MARGIN=1,function(i){ 
#' dat<-rnorm(i)
#' return(c(
#'KL= KL.bmop(object = bmop_fit(dat),dtrue = dnorm),
#'SE= squareError.bmop(object = bmop_fit(dat),dtrue = dnorm)
#' )) })))
#' plot(Ns,res$SE,lwd=3,type="l",log="y",ylab="")
#' points(Ns,res$KL,lwd=3,col="red",type="l")
#' legend("topright",legend = c("Square Error","KL"),
#' col = c("black","red"),lty=1,lwd=3)
#' @docType package
#' @name Rbmop
NULL
#> NULL