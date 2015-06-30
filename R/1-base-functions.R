






#' New bmop object
#'
#' Constructor of bmop object.
#' 
#' @param knots list of numeric vector, knots of the B-spline basis
#' @param order vector order of the B-spline for each variable, values will be
#' recycled eventually
#' @param ctrpoints array of control points
#' @param nk logical, parameter used internally see details on how to use it
#' @return An object of class bmop, a list with elements \code{order}, 
#' \code{knots}, \code{ctrpoints},
#' @details The function create an object of class \code{bmop}, that is, a list
#' with components
#' \itemize{
#' \item \code{order} vector of orders
#' \item \code{knots} list, every element of the list is the vector of knots in
#' one dimension, knots are ordered but can be repeted.
#' \item \code{ctropoints} An array of control points
#' }
#' If \code{nk==FALSE} (default method) the knots provided are sorted and the
#' \code{unique} function is applied. Then the appropriate knots are computed,
#' repeating the first and last knots many times as the order in the considered 
#' dimension, with this choice the bmop object represent a B-spline function
#' regular up to the \code{order-1} derivate in every dimension.
#' If \code{nk==TRUE} the bmop object will be created without any check on the 
#' knots, use carefully this option as it could result in a wrong defined bmop 
#' object.
#' @export
new_bmop<-function(knots,order,ctrpoints = NULL,nk = FALSE){
  d<-length(order)
  l<-length(knots)
  if (l!=d){ order<-rep(order,l)[1:l]
  d<-length(order)
  }
  if (class(knots)=="numeric") { knots<-list(knots)}
  if (!nk){
    knots<-lapply(knots,FUN=unique)
    nknots<-list()
      for (i in 1:d){
        nknots[[i]]<-c(rep(knots[[i]][1],order[i]-1),knots[[i]],
                       rep(knots[[i]][length(knots[[i]])],order[i]-1))
      }

  }
  else {
     nknots<-knots
  }
  if (is.null(ctrpoints)) {
     ctrpoints<-array(0,dim =as.vector(sapply(nknots,FUN=length))-order)
  }
  else{
    ctrpoints<-array(ctrpoints,dim =as.vector(sapply(nknots,FUN=length))-order)
  }
  
  mop<-list()
  mop[["order"]]<-order
  mop[["knots"]]<-nknots
  mop[["ctrpoints"]]<-ctrpoints
  class(mop)<-"bmop"
  return(mop)
}

#' Check if an object's class is bmop
#' 
#' @param object an R object
#' @return logical, TRUE if "bmop" is present among class(object),
#'  FALSE otherwise
#'  @details Just a call of \code{inherits(object,"bmop")}
#' @export
#' @examples
#' bmop<-bmop_fit(rexp(100))
#' is.bmop(bmop)
is.bmop<-function(object){
  inherits(x = object,what = "bmop")
}



#' Evaluation of a bmop object
#'
#' @param x numeric, vector, matrix
#' @param object bmop object
#' @param MIN numeric
#' @return Numeric value (values) of the computed bspline at point
#'  (points) \code{x}.
#'  
#' @details The deBoor algorithm is used, implemented in C.
#' @references Carl de Boor, 
#' On calculating with B-splines, \emph{Journal of Approximation Theory}, 
#' Volume 6, Issue 1, July 1972, Pages 50-62, 
#' \url{http://www.sciencedirect.com/science/article/pii/0021904572900809}.
#' @export
#' @examples
#' bmop<-bmop_fit(rnorm(100))
#' evaluate.bmop(0,bmop)
#' evaluate.bmop(c(-1,0,+1),bmop)
evaluate.bmop<-function(x,object,MIN=bmopPar()$MIN){
   d<-length(object$order)
   if (is.null(dim(x))){
     if (d == 1){
       dim(x)<-length(x)
     }
     if (d > 1){
       dim(x)<-c(1,length(x))
     }
   }
   if (dim(x)[1] != 1) { 
     return(apply(X = x,MARGIN = 1,FUN = evaluate.bmop,object=object,MIN=MIN))
   }
   if (d == 1) {
     dim(x)<-NULL
     return(deboor_c(t = x,k = object$order,knots = object$knots[[1]],
                   ctr = object$ctrpoints,MIN = MIN))
   }
   ctr<-rep(x = 0,times = dim(object$ctrpoints)[1])
   indx<-slice.index(x = object$ctrpoints,MARGIN = 1)
   l<-locate(x[1],object$knots[[1]])
   if (l==0){ return(MIN)}
   if (l>dim(object$ctrpoints)[1]){ l<-dim(object$ctrpoints)[1]}
   for (i in (l-object$order[1] + 1):l){
     m<-new_bmop(knots = object$knots[-1],order = object$order[-1],nk = T)
     o<-object$ctrpoints[indx==i]
     dim(o)=dim(m$ctrpoints)
     m$ctrpoints<-o
     ctr[i]<-evaluate.bmop(object=m,x=x[-1],MIN=MIN)  
   }
   return(deboor_c(t=x[1],k = object$order[1],knots = object$knots[[1]],
                 ctr =ctr ,MIN = MIN))
}



delta<-function(bmop,x,MIN=bmopPar()$MIN){
  d<-length(bmop$order)
  c<-array(dim = dim(bmop$ctrpoints),1)
  for (i in 1:d){
   m<-new_bmop(knots = bmop$knots[i],order = bmop$order[i],nk = T)
   k<-array(dim = dim(bmop$ctrpoints),1)
   ix<-slice.index(x = k,MARGIN = i)
   for (j in 1:(dim(c)[i])){
      m$ctrpoints[j]<-1
      k[ix==j]<-evaluate.bmop(x = x[i],object=m,MIN = MIN)
      m$ctrpoints[j]<-0
  }
  c<-c*k
  }
  return(c)
}


integration_constants<-function(bmop,append.extra=F){
  d<-length(bmop$order)
  dimm<-dim(bmop$ctrpoints)
  if (is.null(dimm)){ dimm<-length(bmop$ctrpoints)}
  c<-array(dim = dimm ,1)
  for (i in 1:d){
    k<-array(dim = dimm,1)
    ix<-slice.index(x = k,MARGIN = i)
    for (j in 1:(dim(c)[i])){
      k[ix==j]<-(bmop$knots[[i]][j+bmop$order[i]]-bmop$knots[[i]][j])/bmop$order[[i]]
    }
    c<-c*k
  } 
  if (append.extra){ bmop$extra$C<-c
              return(bmop)
   }
  return(c)
}
 
#' Plot of bmop object
#' 
#' @param x bmop object
#' @param N positive integer
#' @param type graphical parameter, see \code{\link{par}}
#' @param contour graphical parameter, see \code{\link{par}}
#' @param persp logical
#' @param file optional file name, where to save a pdf copy of the plot 
#' @param MIN parameter of \code{\link{evaluate.bmop}}
#' @param ... graphical parameters as \code{col,main,...}, see 
#' \code{\link{plot}} or \code{\link{filled.contour}}
#' @return plot of 1-d or 2-d (contour and perspective) bmop, an empty plot and 
#' a warning message are returned if the bmop object has more than 2 dimensions.
#' @details For the 2-d persp (\code{persp==TRUE}) plot the \code{rgl} package 
#' is needed.
#' In order to obtain a nicer result in the contour plot, try to supply a 
#' \code{col} or \code{color.palette} parameter, like in the example below 
#' (\code{RColorBrewer} package is needed).
#' @seealso \link{points.bmop}
#' @examples
#' \dontrun{
#' bmop<-bmop_fit(data.frame(rnorm(100),rnorm(100)))
#' colFun<-
#' grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,name="YlGnBu"))
#' plot(bmop,color.palette=colFun)
#' }
#' @export
plot.bmop<-function(x,N=1000,type="l",contour=TRUE,persp=FALSE,file=NULL,MIN=0,
                    ...){
  if (length(x$order) == 1){
    t<-seq(from = x$knots[[1]][1], to = max(x$knots[[1]]),
           by = ( max(x$knots[[1]]) - min(x$knots[[1]]) ) / N )
    y<-sapply(t, FUN = evaluate.bmop, object=x, MIN=MIN)
    if (!is.null(file)){
      pdf(file = file)
    }
    plot(t,y,type = type,...)
    if (!is.null(file)){ dev.off()}
  }
  if (length(x$order) == 2){
    N<-floor(sqrt(N))
    t1<-seq(from = x$knots[[1]][1],to =max(x$knots[[1]])
            ,by = (max(x$knots[[1]])-min(x$knots[[1]]))/N )
    t2<-seq(from = x$knots[[2]][1],to =max(x$knots[[2]])
            ,by = (max(x$knots[[2]])-min(x$knots[[2]]))/N )
    grid<-as.matrix(expand.grid(t1,t2))
    z<-evaluate.bmop(object = x,x = grid,MIN = MIN)
    z<-matrix(nrow=N+1,ncol=N+1,z,byrow=F)
    if (persp){
      if (requireNamespace("rgl", quietly = TRUE)){
      rgl::persp3d(x=t1,y=t2,z=z,col="red",...)
      }
    }
    if (contour) {
      if (!is.null(file)) { pdf(file=file) }
      zlim <- range(z, finite=TRUE)
      zlim[1] <- 0
      nlevels <- 20
      levels <- pretty(zlim, nlevels)
      nlevels <- length(levels)
      filled.contour(x=t1,y=t2,z=z,nlevels=nlevels,levels=levels,las=1,...)
      if (!is.null(file)){ dev.off()}
    }
  }
  if (length(x$order)>  2){
    warning("plotting in more than 3 dimensions is not possible
            , an empty plot is returned")
    return(plot.new())
  }
}

#' Summary of a bmop object
#'
#' @param object bmop object
#' @param ... compatibility with \code{\link{summary}}
#' @export
#' @examples 
#' data(trees)
#' bmop<-bmop_fit(data=trees$Height)
#' summary(bmop)
summary.bmop<-function(object,...){
  summ<-list()
  nknots<-sapply(X = lapply(object$knots,unique),FUN  = length)
  orders<-object$order
  mins<-sapply(X=object$knots,FUN = min)
  maxs<-sapply(X=object$knots,FUN = max)
  summ$numinfo<-data.frame(Dim=1:length(orders),
                           Numb.of.Knots=nknots,Order=orders,Min=mins,Max=maxs)
  summ$logLik<-object$logLik
  class(summ)<-"summary.bmop"
  return(summ)
}


#' @export 
print.summary.bmop<-function(x,...){
 cat("B-MoP object \n")
 cat("\n")
 print.data.frame(x$numinfo,row.names = F,right = T)
 cat("\n")
 if (!is.null(x$logLik)){
  cat("logLik: ")
  cat(x$logLik)
  cat("\n")
 }
}

#' Print bmop objects
#'
#'@param x bmop object
#'@param ... see \code{\link{print}}
#'@export 
print.bmop<-function(x,...){
  print(summary(x))
  cat("\n")
  if (length(x$order<3)){
    cat("Knots: \n")
    sapply(lapply(x$knots,unique),function(x){ cat(x) 
                                               cat("\n")})
    cat("\n")
    cat("Control Points: \n")
    cat(x$ctrpoints)
  }
}

#' Plot points from bmop 
#'
#' @param x bmop object
#' @param N number of points to plot
#' @param ... graphical parameters see \code{\link{par}}
#' @details As \code{\link{points}}, this functions provide a way to add the 
#' plot of a bmop to an existing plot.
#' @seealso \code{\link{plot.bmop}}.
#' @export
#' @examples 
#' data<-rnorm(200)
#' bmopE<-bmop_fit(data)
#' bmopS<-search_bmop(data)
#' plot(bmopE)
#' points(bmopS,type="l",col="red")
points.bmop<-function(x,N=100,...){
  t<-seq(from = x$knots[[1]][1],to =max(x$knots[[1]]),
         by = (max(x$knots[[1]])-min(x$knots[[1]]))/N )
  y<-sapply(t,FUN = evaluate.bmop,object=x)
  points(t,y,...)
}



#' Integrate bmop object over the support  
#'
#' @param object bmop object
#' @export
#' @examples 
#' data<-rnorm(200)
#' bmop<-bmop_fit(data)
#' integrate.bmop(bmop)
integrate.bmop<-function(object){
  return( sum(integration_constants(object)*object$ctrpoints))
  #return(integrate(f = function(x){ sapply(x,FUN=evaluate,object=object)},
  #lower = object$knots[[1]][1],
  #upper = object$knots[[1]][length(object$knots[[1]])]))
}

#' Generate sequence of knots for bmop objects
#'
#' @param data data.frame matrix vector
#' @param N positive integer, the number of knots
#' @param method string "uniform" or "quantiles" 
#' @param Min numeric
#' @param Max numeric
#' @export
#' @examples 
#' data<-rnorm(200)
#' bmopE<-bmop_fit(data)
#' bmopS<-search_bmop(data)
#' plot(bmopE)
#' points(bmopS,type="l",col="red")
generate_knots<-function(data=NULL,N=5,method="uniform",Min=NULL,Max=NULL){
    
    D<-length(N)
    if (inherits(data,"histogram")|inherits(data,"bins")){
      data<-as.data.frame(data$mids)
    }
    if (!is.null(data)){
    data<-as.data.frame(data)
    D<-dim(data)[2]
    if (length(N)!=D){ N<-rep(N,D)[1:D] }
    }
    else { method<-"uniform"}
    if (class(Max)=="list"){Max<-rapply(Max,c)}
    if (class(Min)=="list"){Min<-rapply(Min,c)}
    return(lapply(X = 1:D,FUN = function(i){
      M<-Max[i]
      m<-Min[i]
      if (is.null(Max[i])){ M<-max(data[,i])}
      if (is.null(Min[i])){m<-min(data[,i])}
      if (is.infinite(M)){ M<-Max[i]}
      if (is.infinite(m)){ m<-Min[i]}
      if (method=="uniform"){    
        return(seq(from=m,to=M,by=(M-m)/(N[i])))
      } 
      if (method=="quantiles"){
        return(quantile(x = data[ ,i],probs = (0:N[i]) /N[i]))
      }
    }))
}

#' Add a new knot to a bmop object
#'
#'@param bmop a bmop object
#'@param value  the new knot to be added
#'@param MARGIN the dimension where the new knot has to be added 
#'@return a bmop object
#'@export 
add_knots.bmop<-function(bmop,value,MARGIN=1){
  d<-length(bmop$order)
  knots<-bmop$knots[[MARGIN]]
  l<-locate(value,knots)
  if (l==0){ warning("new knot is outside original interval, original bmop 
                     is returned")
             return(bmop)}
  if (l==(length(knots)-bmop$order[MARGIN]+1)){ warning("new knot is outside 
                                                        original interval, 
                                                        original bmop is
                                                        returned")
                                               return(bmop)}
  knots<-sort(c(knots,value))
  nknots<-bmop$knots
  nknots[[MARGIN]]<-knots
  mop1<-new_bmop(knots = nknots,order = bmop$order,nk = T)
  return(mop1)
}


#' @importFrom stats logLik
NULL



#' Log-Likelihood of bmop object
#'
#' @param object a bmop object
#' @param data matrix, data.frame or vector of observation
#' @param ... some methods require additional arguments, 
#' see \code{\link[stats]{AIC}}
#' @seealso AIC.bmop 
#' @export
logLik.bmop<-function(object,data,...){ 
  if (inherits(data,"histogram")|inherits(data,"bins")){
    counts<-data$counts
    data<-data$mids

    N<-sum(counts)
  }
  else{
    data<-as.data.frame(data)
    counts<-rep(1,dim(data)[1])
    N<-dim(data)[1]
    
  }
  ll<-evaluate.bmop(object = object,x=data,MIN=10^(-10))
  ll<-sum(log(ll)*counts)
  return(ll)
}

#' Dimension of bmop 
#'
#' @param x a bmop object
#' @return the number of free parameters of \code{x}
#' @export
dim.bmop<-function(x){
  return(prod(dim(x$ctrpoints)))
}




#' @importFrom stats BIC AIC
NULL



#' Akaike Information Criteria for bmop 
#'
#' @param object a bmop object
#' @param data data.frame, matrix, vector of observations
#' @param corrected logical 
#' @param ... see \code{\link[stats]{AIC}}
#' @param k penalization coefficients, see \code{\link[stats]{AIC}}
#' @export 
AIC.bmop<-function(object,data,corrected=F,...,k = 2){
  a<-0
  d<-dim.bmop(x =object)
  if (inherits(data,"histogram")|inherits(data,"bins")){
    counts<-data$counts
  }
  else{
    counts<-rep(1,dim(data)[1])
  }
  if (corrected){ a<- 2*d*(d+1)/(sum(counts)-d-1) }
  return(-2*logLik(object,data = data)+k*d+a)
}

#' Bayesian Information Criteria for bmop
#'
#' @param object a bmop object
#' @param data data.frame, matrix, vector of observations
#' @param corrected logical
#' @param ... see \code{\link[stats]{BIC}}
#' @export
BIC.bmop<-function(object,data,corrected=F,...){
  if (inherits(data,"histogram")|inherits(data,"bins")){
    N<-sum(data$counts)
  }
  else{
    N<-dim(as.data.frame(data))[1]
  }
  return(AIC.bmop(object=object,data=data,k=log(N),corrected=corrected))
}

#' Lower Limit of bmop
#'
#'@param object a bmop object
#'@return numeric value
#'@export
lower.bmop<-function(object){
  return(sapply(object$knots,min))
}

#' Upper Limit of bmop
#'
#'@param object a bmop object
#'@return numeric value
#'@export
upper.bmop<-function(object){
  return(sapply(object$knots,max))
}


#' Mean value for a bmop density
#'
#' @param x a bmop object
#' @param ... for compatibility with \code{\link{mean}}
#' @return numeric value, the mean of the bmop density
#' @export
#' @examples 
#' bmop<-bmop_fit(rnorm(100))
#' mean(bmop)
mean.bmop<-function(x,...){
  object<-x
  return(
    cubature::adaptIntegrate(
      f = function(x){return(x*evaluate.bmop(x = x,object = object,MIN=0))},
      lowerLimit = lower.bmop(object),upperLimit = upper.bmop(object))$integral)
}


#' Convert an bmop object to a function
#' 
#' @param x an bmop object
#' @param MIN non negative number, MIN value for evaluation of x
#' @param ... compatibility with \code{\link{as.function}}
#' @export 
#' @examples
#' bmop<-bmop_fit(rnorm(200))
#' fun<-as.function(bmop)
#' fun(0)
#' fun(3)
#' plot(fun,from=-2,to=+2)
as.function.bmop<-function(x,MIN=0,...){
  bmop<-x
  return(function(x,...){return(evaluate.bmop(x = x,object = bmop,MIN = MIN ))})
}
