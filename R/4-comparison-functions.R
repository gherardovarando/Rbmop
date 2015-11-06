
#' Square Error between bmop and true density
#'
#' @param object bmop object
#' @param dtrue function
#' @param ... optional arguments to be passed to \code{dtrue}
#' @return Numeric, the square error between bmop denisty and true density
#' @export
#' @examples 
#' data<-rnorm(200)
#' bmop1<-bmop_fit(data)
#' bmopPar(mle=TRUE)
#' bmop2<-bmop_fit(data)
#' squareError.bmop(bmop1,dtrue=dnorm)
squareError.bmop<-function(object,dtrue,lowerLimit=NULL,upperLimit=NULL,...){
  if (!requireNamespace("cubature", quietly = TRUE)){
    warning("cubature package is required to compare bmop object")
    return(NULL)
  }
  if (is.null(lowerLimit)){
    lowerLimit=lower.bmop(object)
  }
  if (is.null(upperLimit)){
    upperLimit=upper.bmop(object)
  }
  f<-as.function(object)
  ff<-function(x){ (f(x)-dtrue(x,...))^2}
  result<-cubature::adaptIntegrate(f = ff,lowerLimit = lowerLimit,
                         upperLimit = upperLimit)$integral
  return(result)
} 


#' Kullback-Leibler (KL) divergence between bmop and true density
#'
#' @param object bmop object
#' @param dtrue function
#' @param ... optional arguments to be passed to \code{dtrue}
#' @return Numeric, the KL divergence.
#' @export
#' @examples 
#' data<-rnorm(200)
#' bmop1<-bmop_fit(data)
#' bmopPar(mle=TRUE)
#' bmop2<-bmop_fit(data)
#' KL.bmop(bmop1,dtrue=dnorm)
#' KL.bmop(bmop2,dtrue=dnorm)
KL.bmop<-function(object,dtrue,...){
  if (!requireNamespace("cubature", quietly = TRUE)){
    warning("cubature package is required to compare bmop object")
    return(NULL)
  }
  f<-as.function(object)
  ff<-function(x){ 
    d<-dtrue(x,...)
    o<-f(x)
    return(o*log(o/d))
  }
  result<-cubature::adaptIntegrate(f = ff,lowerLimit = lower.bmop(object),
                                   upperLimit = upper.bmop(object))$integral
  return(result)
} 





#' Plot several bmops and true density
#'
#' @param bmop.list list of bmop objects (1 dimensional)
#' @param dtrue function, true density function 
#' @param colors NULL or vector of colors
#' @param lwd graphical parameter see \code{\link{par}}
#' @param type graphical parameter see \code{\link{par}}
#' @param type.true as \code{type}, but just for true density plot
#' @param col.true as \code{col} graphical parameter but for dtrue density plot
#' @param names.bmop vector for the names showed in legend
#' @param legend.display logical, to show or not the legend
#' @param legend.pos string for the position of the legend 
#' \code{"top"}, \code{"left"}, 
#' \code{"topright"}, etc.
#' @param file an optional file name to save the plot
#' @param ylim parameter to set the plot correctly, set to NULL 
#' for automatic setting
#' @param ... more parameters to be passed to \code{dtrue}
#' @return invisible() 
#' @export
#' @examples 
#' data<-rnorm(50)
#' bmop1<-bmop_fit(data)
#' bmopPar(mle=TRUE)
#' bmop2<-bmop_fit(data)
#' bmopPar(mle=FALSE)
#' comparison_plot(list(bmop1,bmop2),dtrue=dnorm,
#'                  names.bmop=c("Fast","MLE"))
comparison_plot<-function(bmop.list,dtrue=NULL,colors=NULL,lwd=3,type="l",
                          type.true="l",col.true="red",
                          names.bmop=1:length(bmop.list),legend.display=T,
                          legend.pos="topleft",file=NULL,ylim=c(0,1),...){
  if (class(bmop.list)=="bmop"){ bmop.list<-list(bmop.list)}
  D<-sapply(X = bmop.list,function(bmop){
    return(length(bmop$order))
  })
  if (any(D>1)){ warning("compare_plot is only available for one dimansional 
                         bmops, an empty plot is returned") 
                 return(plot.new())
  }
  if (is.null(colors)){ colors<-rainbow(length(bmop.list),start=0.25)}
  if (length(colors)!=length(bmop.list)){
    colors<-rainbow(length(bmop.list),start=0.25) }
  min<-sort(sapply(X = bmop.list,function(bmop){
    return(min(bmop$knots[[1]]))
  }),index.return=T)
  max<-sort(sapply(X = bmop.list,function(bmop){
    return(max(bmop$knots[[1]]))
  }),index.return=T,decreasing=T)
  min.value<-min$x[1]
  min.pos<-min$ix[1]
  max.value<-max$x[1]
  max.pos<-max$ix[1]
  if (!is.null(file)){ pdf(file)}
  if (!is.null(dtrue)){
    tt<-seq(from = min.value,to = max.value,by = (max.value-min.value)/1000)
    yy<-dtrue(tt,...)
    plot(tt,yy,type = type.true,lwd=lwd,col = col.true,xlab="x",ylab=""
         ,ylim=ylim)
    points.bmop(x = bmop.list[[1]],col=colors[1],type=type,lwd=lwd)
  }
  else{ 
    tt<-seq(from = min.value,to = max.value,by = (max.value-min.value)/1000)
    plot(tt,rep(0,length(tt)),type = "l",xlab="x",ylab="",ylim = ylim)
  }
if (length(bmop.list)>=1){
 for (i in 1:length(bmop.list)){
   points.bmop(x = bmop.list[[i]],col=colors[i],type=type,lwd=lwd)
 } }
 if (legend.display){
   if (!is.null(dtrue)){
 legend(legend.pos,legend =c("dtrue",as.character(names.bmop)),
        col=c("red",colors),lty=1,lwd=lwd ,cex=0.8)
}
else{
  legend(legend.pos,legend =c(as.character(names.bmop)),
         col=c(colors),lty=1,lwd=lwd ,cex=0.8)
}
}
if (!is.null(file)){ dev.off()}
invisible()
}

#' Envelop of estimated bmop
#' 
#' This function plots various bmop estimations, from differents datasets.
#'
#' @param n number of bmop densities to learn from different datasets
#' @param N number of observations in every dataset
#' @param rtrue function to generate samples, see \code{\link{rnorm}}
#' @param dtrue true density function, see \code{\link{dnorm}}
#' @param fun a learning function of bmop as \code{\link{bmop_fit}}
#' @param lwd graphical par
#' @param type graphical par
#' @param col.true the color for the true density
#' @param ... additional parameters to be passed to \code{rtrue} and
#'  \code{dtrue}
#' @return \code{invisible() }
#' @export
#' @examples  
#' envelope_plot(n=50,N=50,rtrue=rexp,dtrue=dexp)
envelope_plot<-function(n=100,N=50,rtrue=rnorm,fun=bmop_fit,
                        dtrue=dnorm,lwd=3,type="l",col.true="red",...){
  data.list<-lapply(rep(N,n),FUN = rtrue,...)
  bmop.list<-lapply(data.list,FUN=fun)
  min<-sort(sapply(X = bmop.list,function(bmop){
    return(min(bmop$knots[[1]]))
  }),index.return=T)
  max<-sort(sapply(X = bmop.list,function(bmop){
    return(max(bmop$knots[[1]]))
  }),index.return=T,decreasing=T)
  min.value<-min$x[1]
  min.pos<-min$ix[1]
  max.value<-max$x[1]
  max.pos<-max$ix[1]
  Squareerrors<-sapply(bmop.list,squareError.bmop,dtrue,...)
  
  meanSquareError<-mean(Squareerrors)
  sdSquareError<- sd(Squareerrors)
  
  if (!is.null(dtrue)){
    tt<-seq(from = min.value,to = max.value,by = (max.value-min.value)/1000)
    yy<-dtrue(tt,...)
    plot(tt,yy,type = type,lwd=lwd,col = col.true,xlab="x",ylab="",main=
           paste("Smpl (n)",n,
                 ", Obs (N)",N,
             ", MSE",signif(meanSquareError,digits=4),", SD",
                 signif(sdSquareError,digits=4),sep=" "))
    points.bmop(x = bmop.list[[1]],col="grey",type=type,lwd=lwd)
  }
  else{ plot.bmop(x = bmop.list[[1]],col="grey",type=type,lwd=lwd,main=
                    paste("Smpl (n)",n,
                          ", Obs (N)",N,
                          ", MSE",signif(meanSquareError,digits=4),", SD",
                          signif(sdSquareError,digits=4),sep=" "))
  }
  for (i in 2:length(bmop.list)){
    points.bmop(x = bmop.list[[i]],col="grey",type=type,lwd=lwd)
  } 
  if (!is.null(dtrue)){
    tt<-seq(from = min.value,to = max.value,by = (max.value-min.value)/1000)
    yy<-dtrue(tt,...)
    points.bmop(x = bmop.list[[1]],col="grey",type=type,lwd=lwd)
    points(tt,yy,type = type,lwd=lwd,col = col.true,xlab="x",ylab="")
  }
  else{ points.bmop(x = bmop.list[[1]],col="grey",type=type,lwd=lwd)}
  invisible(list(MSE=meanSquareError,SD=sdSquareError))
}


