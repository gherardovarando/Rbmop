




#' The numbers of occurency of an elemnent into a vector
#'
#' @param x value
#' @param v vector
#' @return numeric value giving the number of times an element appears into 
#' the vector 
multiplicity<-function(x,v){
  return(length(v[x==c(x,v)])-1)
}

#' Position of an element into a sorted vector
#'
#' @param a value
#' @param v vector in ascending order (or it will be ordered)
#' @return positive integer, the position i of the first element 
#' of \eqn{v} (\eqn{v_i}) such that \eqn{v_i<=a<v_{i+1}} 
locate<-function(a,v){
  v<-sort(c(a,v))
  return(max((1:(length(v)))[v==a]-1))
}

#' Metropolis-Hasting Sampler
#'
#' @param N positive integer, the number of observations
#' @param d integer >1, the dimension of the random vector
#' @param densit function, the (multi valued) density function
#' @param h positive integer, the jumping width between recorded observations
#' @param M positive integer, the heating phase length
#' @param xstart d-dimensional vector (or NULL) the starting point of the 
#' Markov chain
#' @param max numeric, truncation parameter
#' @param min numeric, truncation parameter
#' @return matrix of observations
#' @export
#' @examples
#' sample1<-sampler_MH(100,1,dnorm)
#' hist(sample1)
#' plot(bmop_fit(sample1))
#' sample2<-sampler_MH(100,1,dnorm,max=0.5,min=-1)
#' hist(sample2)
sampler_MH<-function(N,d,densit,h=3,M=1000,xstart=NULL,max=+Inf,min=-Inf){
   if (is.null(xstart)){ xstart<-rep(0,d) }
     x<-xstart
     sample<-array(dim=c(M + h * N, d))
     for (i in 1:(M + h * N)){
       sample[i, ]<-x
       ok<-F
       while (!ok){
         xnew<-MASS::mvrnorm(n = 1,mu = x,Sigma = diag(x = rep(1,d)))
         if (all(xnew<=max)&all(xnew>=min)) { ok<-T }
       }
       if ((densit(xnew) / densit(x)) > runif(n = 1)) {
         x<-xnew
       } 
     }
   return(sample[seq(from = M + 1,to = M + h * N,by = h), ])
}


#' Fix Data
#'
#' @param data, vector, array, data.frame, matrix
#' @return if \code{is.null(dim(data))} then the function set 
#' \code{dim(data)<-c(length(data),1)} and return \code{data}.
fix_data<-function(data){
 return(as.data.frame(data))
}




