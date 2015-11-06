
#' Examples of bmop density estimations
#'
#' Various examples to show the capabilities of bmop estimations.
#' @param N positive integer, the number of observations
#' @param m1,m2,m3 location parameters
#' @param lambda mixing coefficient, vector or double
#' @return All the functions return an invisible list contating the 
#' generated dataset, the estimated bmops and the true density function.
#' See example on how, for example, plot the default kernel density estimation
#' on the same dataset.
#' If \code{bmopPar(mle=TRUE)} is called before calling an example function 
#' the \code{bmop_fit} function will be computed with maximum likelihood 
#' estimation.
#' @details 
#' This functions generate datasets of \code{N} observations.
#' The function \code{bmop_fit} is then used to estimate the density.
#' A comparative plot is then returned.  
#'  \describe{
#'  \item{\code{ex_bmop_gaussian2Mixture}}{\code{lambda}-mixture 
#'  of two Gaussian densities with unitary variance
#'  and means \code{m1} and \code{m3}.}
#'  \item{\code{ex_bmop_gaussian3Mixture}}{\code{lambda}-mixture 
#'  of three Gaussian densities with unitary variance
#'  and means \code{m1}, \code{m2} and \code{m3}.}
#'  \item{\code{ex_bmop_gaussianBetaGamma}}{\code{lambda}-mixture 
#'  of a Gaussian density with unitary variance
#'  and mean \code{m1}, a Beta density with \code{shape1=2 shape2=5 ncp=m2}  and 
#'  a Gamma density with \code{shape=9 scale=m3/9}.}
#'  }
#' @name Examples_bmop
#' @examples
#' Ex<-ex_bmop_gaussian2Mixture()
#' points(density(Ex$data),type="l",col="blue",lwd=3)
NULL
#> NULL


#' @rdname Examples_bmop
#' @export
ex_bmop_gaussian2Mixture<-function(N=1000,m1=-3,m2=0,lambda=0.5){
  lambda<-abs(lambda)
  if (lambda>1){lambda<-1/lambda}
  D<-c(rnorm(floor(N*lambda),mean=m1),rnorm(floor(N*(1-lambda)),mean=m2))
  f<-function(x){
    return((dnorm(x,mean=m1)*lambda+dnorm(x,mean=m2)*(1-lambda)))
  }
  bmop<-bmop_fit(D)
  bmopS<-search_bmop(D)
  comparison_plot(list(bmop,bmopS),f,names.bmop = c("bmop_fit","search_bmop"),
                  ylim=NULL)
  return(invisible(list(data=D,true=f,bmop=bmop,bmopS=bmopS)))
}


#' @rdname Examples_bmop
#' @export
ex_bmop_gaussian3Mixture<-function(N=1000,m1=-3,m2=0,m3=+3,lambda=c(1,1,1)){
  lambda<-rep(lambda,3)[1:3]
  lambda<-abs(lambda)/sum(lambda)
  D<-c(rnorm(floor(N*lambda[1]),mean=m1),rnorm(floor(N*(lambda[2])),mean=m2),
       rnorm(floor(N*(lambda[3])),mean=m3))
  f<-function(x){
    return((dnorm(x,mean=m1)*(lambda[1])+dnorm(x,mean=m2)*(lambda[2])+
              dnorm(x,mean=m3)*(lambda[3])))
  }
  bmop<-bmop_fit(D)
  bmopS<-search_bmop(D)
  comparison_plot(list(bmop,bmopS),f,names.bmop = c("bmop_fit","search_bmop"),
                  ylim=NULL)
  return(invisible(list(data=D,true=f,bmop=bmop,bmopS=bmopS)))
}


#' @rdname Examples_bmop
#' @export
ex_bmop_gaussianBetaGamma<-function(N=1000,m1=-3,m2=2,m3=3,lambda=c(1,1,1)){
  lambda<-rep(lambda,3)[1:3]
  lambda<-abs(lambda)/sum(lambda)
  D<-c(rnorm(floor(N*lambda[1]),mean=m1),rbeta(floor(N*(lambda[2])),shape1=2,
                                               shape2=5,ncp=m2),
       rgamma(floor(N*(lambda[3])),,shape=3,scale = m3/3))
  f<-function(x){
    return((dnorm(x,mean=m1)*(lambda[1])+
              dbeta(x,,shape1=2,shape2=5,ncp=m2)*(lambda[2])+
              dgamma(x,shape=3,scale = m3/3)*(lambda[3])))
  }

  bmop<-bmop_fit(D)
  bmopS<-search_bmop(D)
  comparison_plot(list(bmop,bmopS),f,names.bmop = c("bmop_fit","search_bmop"),
                  ylim=NULL)
  return(invisible(list(data=D,true=f,bmop=bmop,bmopS=bmopS)))
  
}
