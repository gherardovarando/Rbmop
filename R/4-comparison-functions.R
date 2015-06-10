
#' Compare bmop with true density
#'
#' @param object bmop object
#' @param dtrue function
#' @param measure string, \code{"MSE"}, \code{"MAE"}, \code{"MAX"}
#' @param method string, \code{"grid"} or \code{"montecarlo"}
#' @param densit function or string \code{"uniform"}
#' @param N positive integer, number of points of comparison
#' @param ... optional arguments to be passed to \code{dtrue}
#' @export
#' @examples 
#' data<-rnorm(200)
#' bmop1<-bmop_fit(data)
#' bmopPar(mle=TRUE)
#' bmop2<-bmop_fit(data)
#' compare.bmop(bmop1,dtrue=dnorm)
#' compare.bmop(bmop2,dtrue=dnorm,method="montecarlo")
compare.bmop<-function(object,dtrue,measure="MSE",method="grid",
                       densit=dtrue,N=100,...){
  
  if (method=="grid"){
    data<-expand.grid(lapply(object$knots,
                             FUN = function(x){ return(min(x)+(0:N)*max(x)/N) 
                                                }))
    data<-as.matrix(data)
  }
  if (method=="montecarlo"){
    if (is.character(densit)){ 
      if (densit=="uniform"){ 
        densit<-function(x){ 
          if (all(x<sapply(object$knots,FUN = max)) && 
                all(x>sapply(object$knots,FUN=min)) ){
          return(1) } 
          else { return(0)} 
        } 
      } 
    }
    data<-sampler_MH(N = N^length(object$order),
                     d = length(object$order),densit = densit,h = 3,M = 1000)
  }
  result<-list()
  for (m in measure){
  if (m=="MSE"){ 
    result$MSE<-(sum((evaluate.bmop(object = object,
                              x = data)-dtrue(data,...))^2)/
             (N^length(object$order))) }
  if (m=="MAE"){ 
    result$MAE<-(sum(abs(evaluate.bmop(object = object,x = data)-dtrue(data,...)))
           /(N^length(object$order))) }
  if (m=="MAX"){ 
    result$MAX<-(max(abs(evaluate.bmop(object = object,x = data)-dtrue(data,...)))) }
  }
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
#' @param ... more parameters to be passed to \code{dtrue}
#' @return invisible() 
#' @export
#' @examples 
#' data<-rnorm(200)
#' bmop1<-bmop_fit(data)
#' bmopPar(mle=TRUE)
#' bmop2<-bmop_fit(data)
#' comparison_plot(list(bmop1,bmop2),true=dnorm,
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
#' @return invisible() 
#' @export
#' @examples  
#' envelope_plot(n=50,N=50,rtrue=rexp,dtrue=dexp)
envelope_plot<-function(n=100,N=50,rtrue=rnorm,fun=bmop_fit,
                        dtrue=dnorm,lwd=3,type="l",col.true="red",...){
  data.list<-lapply(rep(N,n),FUN = rtrue)
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
  if (!is.null(dtrue)){
    tt<-seq(from = min.value,to = max.value,by = (max.value-min.value)/1000)
    yy<-dtrue(tt)
    plot(tt,yy,type = type,lwd=lwd,col = col.true,xlab="x",ylab="")
    points.bmop(x = bmop.list[[1]],col="grey",type=type,lwd=lwd)
  }
  else{ plot.bmop(x = bmop.list[[1]],col="grey",type=type,lwd=lwd)
  }
  for (i in 2:length(bmop.list)){
    points.bmop(x = bmop.list[[i]],col="grey",type=type,lwd=lwd)
  } 
  if (!is.null(dtrue)){
    tt<-seq(from = min.value,to = max.value,by = (max.value-min.value)/1000)
    yy<-dtrue(tt)
    points.bmop(x = bmop.list[[1]],col="grey",type=type,lwd=lwd)
    points(tt,yy,type = type,lwd=lwd,col = col.true,xlab="x",ylab="")
  }
  else{ points.bmop(x = bmop.list[[1]],col="grey",type=type,lwd=lwd)}
  invisible()
}


