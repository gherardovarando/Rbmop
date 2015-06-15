





multiplicity<-function(x,v){
  return(length(v[x==c(x,v)])-1)
}


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



fix_data<-function(data){
 return(as.data.frame(data))
}


#' as.bins
#'@export
as.bins<-function(data,breaks=nclass.FD,...){
  bins<-list()
  data<-as.data.frame(data)
  if (dim(data)[2]==1){
    H<-hist(data,breaks=breaks)
    bins$mids<-H$mids
    bins$counts<-H$counts
    class(bins)<-"bins"
    return(bins)
  }
  if (is.function(breaks)){
    Ns<-sapply(data,FUN = breaks)
  }
  else{
    Ns<-rep(breaks,times = dim(data)[2])[1:(dim(data)[2])]
    }
  Seqs<-lapply(1:(dim(data)[2]),function(i){
    return(pretty(x = data[,i],n = Ns[i]))
  })
  Mids<-lapply(Seqs,function(seq){
    mids<-c()
    for (i in 2:length(seq)){
    mids<-c(mids,mean(seq[(i-1):i]))
    }
    return(mids)
  })
  bins$mids<-expand.grid(Mids)
  bins$counts<-c(rep(0,dim(bins$mids)[1]))
  for (i in 1:(dim(data)[1])){
    pos<-sapply(1:(dim(data)[2]),FUN = function(j){
      locate(data[i,j],Seqs[[j]])
    })
    bins$counts[sum((pos)*(Ns)^(0:(length(pos)-1)))]<-
      bins$counts[sum((pos)*(Ns)^(0:(length(pos)-1)))]+1
  }
  idx<-bins$counts!=0
  bins$counts<-bins$counts[idx]
  bins$mids<-bins$mids[idx,]
  class(bins)<-"bins"
  return(bins)
  
}



