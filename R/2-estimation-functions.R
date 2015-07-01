.parnames<-c("mle","N","order","alpha","knotsMethod","k","toll","repMax","MIN"
             ,"autoReduce")
.bmopenv<-new.env()
lockBinding(".bmopenv",environment())
lockBinding(".parnames",environment())
assign(x = ".bmopPars",value = list(mle = FALSE, 
                                    N=NA,
                                    order=3,
                                    alpha=3,
                                    knotsMethod="uniform",
                                    k=2,
                                    toll=10^(-10),
                                    repMax=100,
                                    MIN=10^(-10),
                                    autoReduce=200),
       envir = .bmopenv)
lockEnvironment(env = .bmopenv,bindings = ".bmopPars")

 
#' Set Rbmop parameters
#' 
#' Set appropriate global parameters to be used by \code{Rbmop}
#' functions . This is the appropriate way to change those parameters.
#' The use is similar to \code{par()} in package \code{graphics}.
#' 
#'  @param ...  see details.
#'  
#'  @details This function set parameters used for estimation.The call \code{
#'           bmopPar()} just print the values of the paramters.
#'  
#'  
#'           \code{mle=FALSE}: logical. If use maximum-likelihood estimation of 
#'           bmop coefficient, setting \code{mle=TRUE} just force 
#'           \code{repMax=1}.
#'            
#'           \code{N=NA}: If present, the number of knots in every dimensions. 
#'           Vector of positive integer, if needed values will be recycled.
#'           
#'           \code{order=3}: The order of the B-spline in every dimensions, 
#'           vector of positive integer, if needed values will be recycled.
#'           
#'           \code{alpha=3}: The penalization 
#'                         exponent to compute the number of knots. This is the
#'                         default method to compute the number of knots, 
#'                         with the formula: \eqn{ floor(n^(1/alpha))}, where 
#'                         \eqn{n} is the number of observations in the 
#'                         dataset.
#'                         If \code{!is.na(N)} then the number of knots will be 
#'                         set to \eqn{N^d} where \eqn{d} is the number of
#'                          dimensions
#'                         in the dataset (num. of variables).
#'                         
#'          \code{knotsMethod="uniform"}:
#'           \code{"uniform"} or \code{"quantiles"} knots, 
#'           how knots are computed by \code{\link{generate_knots}}. 
#'           
#'           \code{k=2}: Coefficient of \code{AIC} (penalized likelihood), 
#'           positive integer or \code{"BIC"} string. This is used by 
#'           \code{\link{search_bmop}}.
#'           
#'            \code{toll=10^{-10}}: Tollerance for the increment of the likelihood in
#'                         the mle estimation of the coefficient.
#'                         
#'           \code{repMax=100}: maximum number of iteration in the mle 
#'                              estimation of the coefficient  
#'            
#'          \code{MIN=10^{-10}}: This is not a learning parameter but instead 
#'                               define the \code{MIN} parameter in the 
#'                               evaluation of bmop object. Observe that some 
#'                               functions like \code{\link{logLik}} or 
#'                               \code{plot},
#'                               set this parameter independently.
#'                               
#'          \code{autoReduce=200}: This value set the maximum dimension of an 
#'                                 accepted dataset as raw-data, for larger
#'                                 dataset, functions \code{\link{bmop_fit}}
#'                                  and 
#'                                 \code{\link{search_bmop}}
#'                                  will be applied over the 
#'                                 reduced bins (histogram). Setting it to 
#'                                 \code{Inf} disable this features.
#'           @export
bmopPar<-function(...){
  
  l <- list(...)
  
  old <- get(".bmopPars",envir = .bmopenv)
  if (length(l) == 0){
    return(old)
  }
  unlockBinding(sym = ".bmopPars",env = .bmopenv)
  for (a in names(l)){
    if (a %in% .parnames){
    old[[a]] <- l[[a]]
    }
    else{
      warning(paste(a,"is not a Rbmop parameter"))
    }
  }
  if (is.list(l[[1]])){
    for (a in names(l[[1]])){
      old[[a]] <- l[[1]][[a]]  
    }
      
  }
  assign(".bmopPars",old,envir = .bmopenv)
  lockBinding(sym = ".bmopPars",env = .bmopenv)
}


define_bmop<-function(bmop=NULL,data=NULL,Max=NULL,Min=NULL,
                      N=get(".bmopPars",
                            envir = .bmopenv)$N,
                      order=get(".bmopPars",
                                envir = .bmopenv)$order,
                      alpha=get(".bmopPars",
                                envir = .bmopenv)$alpha,
                      method=get(".bmopPars",
                                 envir = .bmopenv)$knotsMethod,...){
  
  if (is.null(alpha)){alpha<-3}
  
  if (is.null(order[1])){order <- 3}
  if (is.null(method)){method <- "uniform"}
  
  if (is.bmop(bmop)){
    return(bmop)
  }
  if (is.null(data)){ 
    return(
      new_bmop(ctrpoints = 1,
               knots = generate_knots(data = data,N = max(1,N),
                                      Max = Max,Min = Min),
               order = order))}
  if (inherits(data,what = "histogram")|inherits(data,what = "bins")){
    counts <- data$counts
    data <- data$mids
    data <- as.data.frame(data)
    alpha <- 1/log(sum(counts) ^ ( 1/alpha ),base = dim(data)[1])
  }
  else if (inherits(data,"values")){
    data <- data$mids
    data <- as.data.frame(data)
  }
  else{
  data<-as.data.frame(data)
  }
  if (is.na(N[1])){
    N<-max(1,floor(dim(data)[1] ^ (1 / (dim(data)[2]*alpha))))
    if (length(order)!=dim(data)[2]){ 
      order<-rep(order,dim(data)[2])[1:(dim(data)[2])]}
  }
  return(
  new_bmop(ctrpoints = 1,knots = 
            generate_knots(data = data,N =N,Max = Max,Min=Min,method = method )
            ,order = order))
  }


#'Estimation of bmop density or conditional density
#'
#'@param data data.frame, matrix, vector or object of class "histogram" or
#' "bins"
#' @param density logic, if \code{TRUE} a density is learned
#' @param conditional logic, if \code{TRUE} a conditional density is learned
#' @param Min vector of lower bounds
#' @param Max vector of upper bounds
#' @param ... see \code{bmopPar}
#' @return a bmop object, a density if \code{conditional=FALSE} or a 
#' conditional density if \code{conditional=TRUE}. In the latter case the first 
#' variable of the dataset is considered asa de conditioned one and the rest of 
#' the variables as the conditioning ones. If the dataset has only one variable 
#' a normal density is generated discarding the value of \code{conditional}. 
#' @export
#' @examples
#' plot(bmop_fit(rnorm(100)))
#' plot(bmop_fit(hist(rnorm(100000))))
#' #############################
#' Data<-data.frame(rnorm(100),rexp(100))
#' bmop<-bmop_fit(Data)
#' plot(bmop)
#' #############################
#' X<-rnorm(100)
#' Y<-rnorm(100,mean=X)
#' Data<-data.frame(X,Y)
#' bmopPar(mle=TRUE)
#' bmopC<-bmop_fit(Data,conditional=TRUE)
#' plot(bmopC)
bmop_fit<-function(data,density=TRUE,conditional=FALSE,Min=NULL,Max=NULL,...){
  UseMethod("bmop_fit",object = data)
}

#'Estimation of bmop density or conditional density
#'
#'@param data data.frame, matrix, vector or object of class "histogram" or
#' "bins"
#' @param density logic, if \code{TRUE} a density is learned
#' @param conditional logic, if \code{TRUE} a conditional density is learned
#' @param Min vector of lower bounds
#' @param Max vector of upper bounds
#' @param ... see \code{bmopPar}
#' @return a bmop object
#' @export
 bmop_fit.default<-function(data,density=TRUE,conditional=FALSE,Min=NULL,
                            Max=NULL,...){
 return(bmop_fit.data.frame(as.data.frame(data),density=density,
                            conditional=conditional,
                     Min=Min,Max=Max,...) )
}


  
  
#'Estimation of bmop density or conditional density
#'
#'@param data data.frame, matrix or vector
#'the variables must be in the right order (the columns of data)
#' @param density logic, if \code{TRUE} a density is learned
#' @param conditional logic, if \code{TRUE} a conditional density is learned
#' @param Min vector of lower bounds
#' @param Max vector of upper bounds
 #' @param bmop  a bmop object
#' @param ... see \code{bmopPar}
#' @return a bmop object
#'@export
bmop_fit.data.frame<-function(data,density=TRUE, conditional=F,
                              Min=NULL,
                              Max=NULL,
                              bmop=NULL,
                              ... ){
  if ((!density)&(!conditional)){
    if (dim(data)[2]==1){
      density=TRUE
    }
    else{
    H<-list()
    class(H)<-"values"
    H$mids<-data[,-(dim(data)[2])]
    H$counts<-data[,(dim(data)[2])]
    return(bmop_fit(H,Min=Min,Max=Max,bmop=bmop,...))
    }
  }
  data<-as.data.frame(data)
  if (dim(data)[1]>bmopPar()$autoReduce){
    warning("Data dimension exceded bmopPar()$autoreduce parameter. 
            The data is grouped into bins. Modifying bmopPar()autoReduce 
            to a greater value (or Inf) prevent that behaviour.
            See the help for more
            informations.")
    if (conditional){
      breaks<-max(1,floor(nclass.FD(data[,1])^{1/(dim(data)[2])}) )
    }
    else{
      breaks<-max(1, floor(max(sapply(data,nclass.FD))^{1/(dim(data)[2])} ))
    }
    return(bmop_fit.bins(data = as.bins(data = data,breaks = breaks),
                                        density=density,
                                        conditional=conditional,...))
  }
  m<-define_bmop(bmop = bmop,data = data,Max = Max,Min = Min,...)
  m<-normalize.bmop(m)
  d<-length(m$order) 
  C<-integration_constants(m)
  N<-dim(data)[1]
  D<-apply(data,MARGIN = 1,FUN = delta,bmop=m,MIN=10^(-10)) 
  dim(D)<-c(dim(m$ctrpoints),dim(data)[1])
  E<-slice.index(D,MARGIN = length(dim(D)))
  m$logLik<-logLik.bmop(object=m,data=data)
  
  if (!bmopPar()$mle){ repmax <- 1}
  else{
    repmax<-bmopPar()$repMax
  }
  
  if (!density){
    repmax <- 1 
    # new bmop_fit.sskskfklak
  }
  
  for (s in 1:repmax){
    K<-array(dim=dim(m$ctrpoints),0)
    
      for (i in 1:(dim(data)[1])){
        K<-K+D[E==i]/sum(D[E==i]*m$ctrpoints)
      }
    
    
    mnew<-m
    mnew$ctrpoints<-m$ctrpoints*K
    KK<-N
    if (density){
      KK<- N*C
    }
    if (conditional){
      C<-integration_constants(new_bmop(knots = m$knots[1],
                                        order = m$order[1]))
      ix<-slice.index(mnew$ctrpoints,MARGIN=1)
      KK <- mnew$ctrpoints
      for (i in 1:(dim(mnew$ctrpoints)[1])){
        KK [ix==i] <- (C[i]*apply(mnew$ctrpoints,MARGIN = -1,sum)) 
      }
    }
    mnew$ctrpoints<- mnew$ctrpoints/ KK
    mnew$logLik<-logLik.bmop(object = mnew,data = data)
    if ((mnew$logLik-m$logLik)<bmopPar()$toll){
      break 
    }
    m<-mnew
  }
  
  return(mnew)
}


#'Estimation of bmop density or conditional density
#'
#'@param data \code{histogram} or \code{bins} object
#'the variables must be in the right order
#' @param density logic, if \code{TRUE} a density is learned
#' @param conditional logic, if \code{TRUE} a conditional density is learned
#' @param Min vector of lower bounds
#' @param Max vector of upper bounds
#' @param bmop  a bmop object
#' @param ... see \code{bmopPar}
#' @return a bmop object
#'@export
bmop_fit.bins<-function(data, density=TRUE, conditional=F,
                              Min=NULL,
                              Max=NULL,
                              bmop=NULL,
                              ... ){
  if ((!density)&(!conditional)){
   class(data)<-"values"
   return(bmop_fit(data,Min=Min,Max=Max,bmop=bmop,...))
  }
  m<-define_bmop(bmop = bmop,data = data,Max = Max,Min = Min,...)
  m<-normalize.bmop(m)
  d<-length(m$order) 
  C<-integration_constants(m)
  counts<-data$counts
  DD<-data
  data<-as.data.frame(data$mids)
  N<-sum(counts)
  D<-apply(data,MARGIN = 1,FUN = delta,bmop=m,MIN=10^(-10)) 
  dim(D)<-c(dim(m$ctrpoints),dim(data)[1])
  E<-slice.index(D,MARGIN = length(dim(D)))
  m$logLik<-logLik.bmop(object=m,data=data)
  
  if (!bmopPar()$mle){ repmax <- 1}
  else{
    repmax<-bmopPar()$repMax
  }
  for (s in 1:repmax){
    K<-array(dim=dim(m$ctrpoints),0)
    for (i in 1:(dim(data)[1])){
     K<-K+D[E==i]*counts[i]/sum(D[E==i]*m$ctrpoints)
     # K<-K+D[E==i]*counts[i]
    }
    mnew<-m
    mnew$ctrpoints<-m$ctrpoints*K
    KK<-dim(data)[1]*C
    if (density){
      KK<- N*C
    }
    if (conditional){
      C<-integration_constants(new_bmop(knots = m$knots[1],
                                        order = m$order[1]))
      ix<-slice.index(mnew$ctrpoints,MARGIN=1)
      KK <- mnew$ctrpoints
      for (i in 1:(dim(mnew$ctrpoints)[1])){
        KK [ix==i] <- (C[i]*apply(mnew$ctrpoints,MARGIN = -1,sum)) 
      }
    }
    mnew$ctrpoints<- mnew$ctrpoints/ KK
    mnew$logLik<-logLik.bmop(object = mnew,data = DD)
    if ((mnew$logLik-m$logLik)<bmopPar()$toll){
      break 
    }
    m<-mnew
  }
  mnew$info$density<- density
  mnew$info$conditional <- conditional
  return(mnew)
}

#'Estimation of bmop density or conditional density
#'
#'@param data \code{histogram} or \code{bins} object
#'the variables must be in the right order
#' @param density logic, if \code{TRUE} a density is learned
#' @param conditional logic, if \code{TRUE} a conditional density is learned
#' @param Min vector of lower bounds
#' @param Max vector of upper bounds
#' @param bmop  a bmop object
#' @param ... see \code{bmopPar}
#' @return a bmop object
#'@export
bmop_fit.histogram <- bmop_fit.bins

#'Estimation of bmop regression fucntions
#'
#'@param data \code{values} object
#'the variables must be in the right order
#' @param Min vector of lower bounds
#' @param Max vector of upper bounds
#' @param bmop  a bmop object
#' @param ... see \code{bmopPar}
#' @return a bmop object
#'@export
bmop_fit.values<-function(data,
                        Min=NULL,
                        Max=NULL,
                        bmop=NULL,
                        ... ){
  m<-define_bmop(bmop = bmop,data = data,Max = Max,Min = Min,...)
 # m<-normalize.bmop(m)
  d<-length(m$order) 
  C<-integration_constants(m)
  counts<-data$counts
  data<-as.data.frame(data$mids)
  D<-apply(data,MARGIN = 1,FUN = delta,bmop=m,MIN=10^(-10)) 
  dim(D)<-c(dim(m$ctrpoints),dim(data)[1])
  E<-slice.index(D,MARGIN = length(dim(D)))
 
    K<-array(dim=dim(m$ctrpoints),0)
    for (i in 1:(dim(data)[1])){
      #K<-K+D[E==i]*counts[i]/sum(D[E==i]*m$ctrpoints)
       print(D[E==i])
       K<-K+D[E==i]*counts[i]
    }
    mnew<-m
    mnew$ctrpoints<- K
    KK<-(dim(data)[1])*C
    mnew$ctrpoints<- mnew$ctrpoints/ KK
  print(C)
  return(mnew)
}

#' Greedy penalized log-likelihood search
#'
#'Aproximation of a density \eqn{f(x_1,\ldots,x_n)} or conditional density
#'@param data data.frame, matrix or vector, the variables must be in
#' the right order (the columns of data)
#' @param conditional logical 
#'@param k positive number or \code{"BIC"}
#'@param corrected logical
#'@param knotsMethod the method to use in knots generation
#'@param ... additional parameters 
#'@return A bmop object, the aproximations of \eqn{f}.
#'@export
#'@examples
#' data<-rnorm(100)
#' bmopS<-search_bmop(data=data)
#' plot(bmopS)
search_bmop<-function(data,conditional=F,k=Rbmop::bmopPar()$k,corrected=FALSE,
                      knotsMethod=Rbmop::bmopPar()$knotsMethod,...){
  
  oldpar<-bmopPar()
  if (inherits(data,"histogram")|inherits(data,"bins")){
    N<-sum(data$counts)
    d<-dim(as.data.frame(data$mids))[2]
    Ddata<-as.data.frame(data$mids)
  }
  else {
    data<-as.data.frame(data)
    Ddata<-data
    N<-dim(data)[1]
    d<-dim(data)[2]
  }
  
  N<-rep(2,d)
  order=rep(3,d)
  bmopPar(N=N,order=order,mle=T)
  mop<-bmop_fit(data = data,conditional = conditional)
  score<-AIC.bmop(object = mop,k = k,data = data,corrected = corrected)
  finish<-F
  while (!finish){
    Mlist<-list()
    finish<-T
    for (i in 1:d){
      norder<-order
      norder[i]<-norder[i]+1
      bmopPar(order=norder,N=N)
      Mlist[[i]]<-bmop_fit(data = data,
                           conditional = conditional)
      Nnew<-N
      Nnew[i]<-Nnew[i]+1
      bmopPar(order=order,N=Nnew)
      Mlist[[d+i]]<-bmop_fit(data=data,
                             conditional=conditional)
    }
    Mlist<-Mlist[!sapply(Mlist,is.null)]
    nScores<-sort(x = sapply(X =Mlist ,FUN = AIC.bmop,data=data,k=k,
                             corrected=corrected),index.return=T)
    if (nScores$x[1]<score){
      mop<-Mlist[[nScores$ix[1]]]
      score<-nScores$x[1]
      N<-sapply(mop$knots,function(x){length(unique(x))})
      order<-mop$order
      finish<-F
    }
  }
  mop$extra$AIC$value<-score
  mop$extra$AIC$k<-k
  bmopPar(oldpar)
  return(mop)
}


