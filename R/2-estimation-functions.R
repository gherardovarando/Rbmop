.parnames<-c("ml","N","order","alpha","knotsMethod","k","MIN")
.bmopenv<-new.env()
lockBinding(".bmopenv",environment())
lockBinding(".parnames",environment())
assign(x = ".bmopPars",value = list(mle = FALSE, 
                                    N=NA,
                                    order=3,
                                    alpha=3,
                                    knotsMethod="uniform",
                                    k=2,
                                    MIN=10^(-10) ),
       envir = .bmopenv)
lockEnvironment(env = .bmopenv,bindings = ".bmopPars")

 
#' Set fitting parameters
#' 
#' Set appropriate parameters to be used by fitting
#' functions.
#' 
#'  @param ... see details.
#'  @details This function set parameters used for estimation.
#'           \code{N}:  The number of knots in every dimensions.
#'           \code{order}: The order of the B-spline in every dimensions.
#'           \code{alpha}: The exponent to compute the number of knots.
#'          \code{knotsMethod}: \code{"uniform"} or \code{"quantiles"} knots. 
#'           
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
  assign(".bmopPars",old,envir = .bmopenv)
  lockBinding(sym = ".bmopPars",env = .bmopenv)
}

#' define bmop, internal function
#'
define_bmop<-function(bmop=NULL,data=NULL,Max=NULL,Min=NULL,
                      N=get(".bmopPars",
                            envir = .bmopenv)$N,
                      order=get(".bmopPars",
                                envir = .bmopenv)$order,
                      alpha=get(".bmopPars",
                                envir = .bmopenv)$alpha,
                      method=get(".bmopPars",
                                 envir = .bmopenv)$knotsMethod){
  
  if (is.null(alpha)){alpha<-3}
  
  if (is.null(order)){order <- 3}
  if (is.null(method)){method <- "uniform"}
  
  if (is.bmop(bmop)){
    return(bmop)
  }
  if (is.null(data)){ 
    return(
      new_bmop(knots = generate_knots(data = data,N = N,Max = Max,Min = Min),
               order = order))}
  data<-fix_data(data)
  if (is.na(N)){
    N<-floor(dim(data)[1] ^ (1 / (dim(data)[2]*alpha)))
    if (length(order)!=dim(data)[2]){ 
      order<-rep(order,dim(data)[2])[1:(dim(data)[2])]}
  }
  return(
  new_bmop(knots = 
            generate_knots(data = data,N =N,Max = Max,Min=Min,method = method )
          ,order = order))
  }



#' Fast Estimation of bmop density
#'
#'Aproximation of a density \eqn{f(x_1,\ldots,x_n)} 
#'@param data data.frame, matrix or vector, 
#'the variables must be in the right order (the columns of data)
#'@param Min real or vector  
#'@param Max real or vector
#'@param ... see \code{bmopPar}
#'@param D internal parameter
#'@return A bmop object, the aproximations of \eqn{f}.
#'@export
#'@examples
#' require("MASS")
#' data<-mvrnorm(n=50,mu=c(0,0),Sigma=diag(x=c(1,1),nrow = 2,ncol = 2))
#' bmop<-estimate_bmop(data=data)
#' plot(bmop,N=20,persp=TRUE)
#' 
#' data<-rnorm(n=1000)
#' bmop<-estimate_bmop(data=data)
#' plot(bmop)
#' summary(bmop)
#' compare.bmop(object=bmop,dtrue=rnorm,N=1000)
estimate_bmop<-function(data,bmop=NULL,conditional=F,Min=NULL,Max=NULL,...,D=NULL){
  mop<-define_bmop(bmop = bmop,data = data,Max = Max,Min=Min,...)
  data<-fix_data(data)
  C<-integration_constants(mop)
  if (is.null(D)){
    D<-apply(data,MARGIN = 1,FUN = delta,bmop=mop) 
    dim(D)<-c(dim(mop$ctrpoints),dim(data)[1])
  }
  E<-slice.index(D,MARGIN = length(dim(D)))
  K<-array(dim=dim(mop$ctrpoints),0)
  for (i in 1:dim(data)[1]){
    K<-K+D[E==i]
  }
  if (conditional){
    C<-integration_constants(new_bmop(knots = mop$knots[1],
                                      order = mop$order[1]))
    ix<-slice.index(mop$ctrpoints,MARGIN=1)
    M<-K
    for (i in 1:(dim(mop$ctrpoints)[1])){
      M[ix==i]<-C[i]*apply(K,MARGIN = -1,sum)
    }
    
  }
  else { M<-C*(dim(data)[1])}
  mop$ctrpoints<-K/M
  mop$info$conditional<-conditional
  return(mop)  
}


#' Fast Estimation of bmop conditional density
#'
#'Aproximation of a conditional density \eqn{f(x| x_1,\ldots,x_n)} 
#'@param data data.frame, matrix or vector, 
#'the variables must be in the right order (the columns of data)
#'@param Min real or vector  
#'@param Max real or vector
#'@param ... see \code{bmopPar}
#'@return A bmop object, the aproximations of \eqn{f}.
#'@export
#'@examples
#' require("MASS")
#' X<-rnorm(100)
#' Y<-rnorm(100,mean=X)
#' data<-data.frame(Y,X)
#' condbmop<-estimate_bmop(data=data)
#' plot(condbmop,N=20,persp=TRUE)
cond_estimate_bmop<-function(data,bmop=NULL,Min=NULL,Max=NULL,...){
  mop<-define_bmop(bmop = bmop,data = data,Max = Max,Min=Min,...)
  data<-fix_data(data)
  ix<-slice.index(mop$ctrpoints,MARGIN=1)
  C<-integration_constants(new_bmop(knots = mop$knots[1],order = mop$order[1]))
  D<-apply(data,MARGIN = 1,FUN = delta,bmop=mop) 
  dim(D)<-c(dim(mop$ctrpoints),dim(data)[1])
  E<-slice.index(D,MARGIN = length(dim(D)))
  K<-array(dim=dim(mop$ctrpoints),0)
  M<-K
  for (i in 1:(dim(data)[1])){
    #dens<-sum(D[E==i])
    K<-K+D[E==i]
  }
  for (i in 1:(dim(mop$ctrpoints)[1])){
    M[ix==i]<-C[i]*apply(K,MARGIN = -1,sum)
  }
  K<-K/M
  mop$ctrpoints<-K
  class(mop)<-c("bmop","cond_bmop")
  return(mop) 
}


#' Maximum conditional-likelihood estimator for the coefficient of B-spline
#'
#'Aproximation of a conditional density \eqn{f(x| x_1,\ldots,x_n)} by maximum 
#'cond-likelihood, the type of mop, knots and order is fixed by the user
#'@param data data.frame, matrix or vector first column is the
#' x variable and the following columns are for 
#' conditioning variable \eqn{x_1,...,x_n}
#'@param bmop bmop object
#'@param Min real or vector  
#'@param Max real or vector
#'@param toll numeric tollerance parameter for the loglik increase
#'@param repmax positive integer maximum number of repetitions
#'@param fresh.start logical if TRUE the cotrol points of mop are set to uniform
#'@param ... see \code{bmopPar}
#'@return if return.mop=TRUE a bmop object, the aproximations of \eqn{f},
#' else a data.frame of control points. 
#'@export
cond_mle_fit_bmop<-function(data,toll=10^(-1),repmax=50,fresh.start=TRUE,
                            bmop=NULL,Min=NULL,Max=NULL,return.bmop=TRUE,...){
  
mop<-define_bmop(bmop = bmop,data = data,Max = Max,Min=Min, ...)
  data<-fix_data(data)
  d<-length(mop$order) 
  C<-integration_constants(new_bmop(knots = mop$knots[1],order = mop$order[1]))
  if (fresh.start){
    a<-array(dim=dim(mop$ctrpoints),1)
    a<-a/sum(C)
  }
  else { a<-mop$ctrpoints}
  N<-dim(data)[1]
  D<-apply(X = data,MARGIN = 1,FUN = delta,bmop=mop,MIN=10^(-20)) 
  dim(D)<-c(dim(mop$ctrpoints),N)
  E<-slice.index(D,MARGIN = length(dim(D)))
  ix<-slice.index(a,MARGIN=1)
  m<-mop
  m$ctrpoints<-a
  for (s in 1:repmax){
    K<-array(dim=dim(a),0)
    M<-K
    for (i in 1:N){
      dens<-sum(D[E==i]*a)
      K<-K+D[E==i]/dens
    }
    anew<-a*K
    for (i in 1:(dim(a)[1])){
      M[ix==i]<-C[i]*apply(anew,MARGIN = -1,sum)
    }
    anew<-anew/M
    mnew<-m
    mnew$ctrpoints<-anew
    if (logLik.bmop(object = mnew,data = data)-
          logLik.bmop(object=m,data=data)<toll){
      if (return.bmop) {
        class(mnew)<-c("bmop","cond_bmop")
        return(mnew)}
      return(anew)
    }
    a<-anew
    m<-mnew
  }
  if (return.bmop){
    class(mnew)<-c("bmop","cond_bmop")
    return(mnew)}
  return(anew)
}


#' Maximum likelihood estimator for the coefficient of B-spline
#'
#'Aproximation of a density \eqn{f(x_1,\ldots,x_n)} by maximum likelihood, 
#'the type of mop, knots and order is fixed by the user
#'@param data data.frame, matrix or vector, 
#'the variable must be in the right order (the columns of data)
#'@param bmop bmop object (optional)
#'@param toll numeric tollerance parameter for the loglik increase
#'@param repmax positive integer maximum number of repetitions
#'@param fresh.start logical if \code{TRUE} the cotrol points of mop are set 
#'                           to uniform
#'@param return.bmop logical, if \code{TRUE} the function return a bmop object
#'@param ... see \code{bmopPar}
#'@return if \code{return.bmop=TRUE} a bmop object,
#' the aproximations of \eqn{f}, else a data.frame of control points. 
#'@export
#'@examples
#' require("MASS")
#' data<-mvrnorm(n=50,mu=c(0,0),Sigma=diag(x=c(1,1),nrow = 2,ncol = 2))
#' bmop<-new_bmop(knots=generate_knots(data),order=c(3,3))
#' bmop<-mle_fit_bmop(data=data,bmop=bmop,return.bmop=TRUE)
#' plot(bmop,N=20,persp=TRUE)
mle_fit_bmop<-function(data,toll=10^(-2),repmax=50,fresh.start=TRUE,
                       return.bmop=TRUE,bmop=NULL,Min=NULL,Max=NULL,...){
  mop<-define_bmop(bmop = bmop,data = data,Max = Max,Min = Min,...)
  data<-fix_data(data)
  d<-length(mop$order) 
  C<-integration_constants(mop)
  m<-mop
  if (fresh.start){
    a<-array(dim=dim(mop$ctrpoints),1)
    a<-a/sum(C)
    m$ctrpoints<-a
  }
  a<-m$ctrpoints
  N<-dim(data)[1]
  D<-apply(data,MARGIN = 1,FUN = delta,bmop=m,MIN=10^(-10)) 
  dim(D)<-c(dim(m$ctrpoints),N)
  E<-slice.index(D,MARGIN = length(dim(D)))
  ll<-logLik.bmop(object=m,data=data)
  for (s in 1:repmax){
    K<-array(dim=dim(a),0)
    h<-(1:N)
    dim(h)<-N
    for (i in 1:N){
      K<-K+D[E==i]/sum(D[E==i]*a)
    }
    anew<-a*K/(N*C)
    mnew<-m
    mnew$ctrpoints<-anew
    llnew<-logLik.bmop(object = mnew,data = data)
    mnew$logLik<-llnew
    if (llnew-ll<toll){
      if (return.bmop) {return(mnew)}
      return(anew)
    }
    a<-anew
    m<-mnew
    ll<-llnew
  }
  if (return.bmop) {return(mnew)}
  return(anew)
}







#' Greedy penalized log-likelihood search
#'
#'Aproximation of a density \eqn{f(x_1,\ldots,x_n)} 
#'@param data data.frame, matrix or vector, the variables must be in
#' the right order (the columns of data)
#'@param k positive number or \code{"BIC"}
#'@param corrected logical
#'@return A bmop object, the aproximations of \eqn{f}.
#'@export
#'@examples
#' require("MASS")
#' data<-rnorm(100)
#' bmopS<-search_bmop(data=data)
#' plot(bmopS)
search_bmop<-function(data,k=Rbmop::bmopPar()$k,corrected=FALSE,
                      knotsMethod=Rbmop::bmopPar()$knotsMethod,...){
data<-fix_data(data)
  N=dim(data)[1]
  if (k=="BIC"){
    k<-log(N)
  }
  d<-dim(data)[2]
  mop<-new_bmop(knots = generate_knots(data = data,N = rep(1,d),
                                       method = knotsMethod),order = rep(2,d))
  mop<-mle_fit_bmop(data = data,bmop = mop)
  score<-AIC.bmop(object = mop,k = k,data = data,corrected = corrected)
  finish<-F
  while (!finish){
    Mlist<-list()
    finish<-T
    for (i in 1:d){
      norder<-mop$order
      norder[i]<-norder[i]+1
      Mlist[[i]]<-new_bmop(knots = mop$knots,order = norder,nk = F)
      Mlist[[i]]<-mle_fit_bmop(data = data,bmop = Mlist[[i]],fresh.start = T)
      Nnew<-sapply(mop$knots,FUN = length)
      Nnew[i]<-Nnew[i]+1
      Mlist[[d+i]]<-new_bmop(knots=generate_knots(data = data,N =Nnew,
                                                  method=knotsMethod),
                                                   order = mop$order)
      Mlist[[d+i]]<-mle_fit_bmop(data=data,bmop=Mlist[[d+i]],fresh.start=T)
    }
    Mlist<-Mlist[!sapply(Mlist,is.null)]
    nScores<-sort(x = sapply(X =Mlist ,FUN = AIC.bmop,data=data,k=k,
                             corrected=corrected),index.return=T)
    if (nScores$x[1]<score){
      mop<-Mlist[[nScores$ix[1]]]
      score<-nScores$x[1]
      finish<-F
    }
  }
  mop$extra$AIC$value<-score
  mop$extra$AIC$k<-k
  return(mop)
}


#' Greedy penalized conditional log-likelihood search
#'
#'Aproximation of a conditional density density \eqn{f(x|x_1,\ldots,x_n)} 
#'@param data data.frame, matrix or vector, the variables 
#'must be in the right order (the columns of data)
#'@param k positive number or \code{"BIC"}
#'@param corrected logical
#'@return A bmop object, the aproximations of \eqn{f}.
#'@export
#'@examples
#' X<-rnorm(100)
#' Y<-rnorm(100,mean=X)
#' data<-data.frame(Y,X)
#' condbmop<-cond_search_bmop(data=data,k="BIC")
#' plot(condbmop,N=20,persp=TRUE)
cond_search_bmop<-function(data,k=Rbmop::bmopPar()$k,corrected=FALSE,
                           knotsMethod=Rbmop::bmopPar()$knotsMethod, ...){
  data<-fix_data(data)
  N=dim(data)[1]
  if (k=="BIC"){
    k<-log(N)
  }
  d<-dim(data)[2]
  mop<-new_bmop(knots = generate_knots(data = data,N = rep(2,d),
                                       method = knotsMethod),order = rep(2,d))
  mop<-cond_mle_fit_bmop(data = data,bmop = mop)
  score<-AIC(object = mop,k = k,data = data,corrected = corrected)
  finish<-F
  while (!finish){
    Mlist<-list()
    finish<-T
    for (i in (1:d)){
      norder<-mop$order
      norder[i]<-norder[i]+1
      Mlist[[i]]<-new_bmop(knots = mop$knots,order = norder,nk = F)
      Mlist[[i]]<-cond_mle_fit_bmop(data = data,bmop = Mlist[[i]],
                                    fresh.start = T)
      Nnew<-sapply(mop$knots,FUN = length)
      Nnew[i]<-Nnew[i]+1
      Mlist[[d+i]]<-new_bmop(knots=generate_knots(data = data,
                                                  N =Nnew,method=knotsMethod),
                                    order = mop$order)
      Mlist[[d+i]]<-cond_mle_fit_bmop(data=data,bmop=Mlist[[d+i]],fresh.start=T)
    }
    Mlist<-Mlist[!sapply(Mlist,is.null)]
    nScores<-sort(x = sapply(X =Mlist ,FUN = AIC.bmop,
                             data=data,k=k,corrected=corrected),index.return=T)
    if (nScores$x[1]<score){
      mop<-Mlist[[nScores$ix[1]]]
      score<-nScores$x[1]
      finish<-F    
    }
  }
  mop$extra$AIC$value<-score
  mop$extra$AIC$k<-k
  class(mop)<-c("bmop","cond_bmop")
  return(mop)
}



