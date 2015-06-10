
#' Clean a bmop object 
#' 
#' remove extra informations from bmop object
#' 
#' @param object a bmop object
#' @return a bmop object
#' @export
clean.bmop<-function(object){
  object$extra<-NULL
  object$logLik<-NULL
  object$AIC<-NULL
  object$info<-NULL
  return(object)
}


#' Marginalize a bmop
#' 
#' marginalize a bmop object, that is integrate a bmop multivariate density 
#' over some of its dimension.
#' @param object a bmop object
#' @param MARGIN  positive integer or vector of positive integer, 
#' the dimensions that has to be marginalized 
#' @param ... additional parameters
#' @return a bmop object over a space of dimension \code{length(MARGIN)} 
#' the result of integrating over the \code{-MARGIN}
#' @export
#' @examples
#' data<-data.frame(rnorm(100),rnorm(100))
#' bmop2d<-bmop_fit(data)
#' bmop1d<-bmop_fit(data[,1])
#' bmop1dmargin<-marginalize.bmop(bmop2d,MARGIN=1)
#' comparison_plot(list(bmop1d,bmop1dmargin),true=dnorm,
#' names.bmop=c("direct est.","marginalized"))
marginalize.bmop<-function(object,MARGIN=1,...){
  mop<-new_bmop(knots = object$knots[MARGIN],order = object$order[MARGIN])
  C<-integration_constants(bmop = object)
  c<-integration_constants(bmop=mop)
  mop$ctrpoints<-apply(C*object$ctrpoints,MARGIN = MARGIN,FUN = sum)
  mop$ctrpoints<-mop$ctrpoints/c
  return(mop)
}


#' Normalize a bmop
#' 
#' Normalize a bmop object, makes it integrate to one.
#' @param object a bmop object
#' @param ... additional parameters
#' @return a bmop object proportional to object, 
#' but such that integrates to one.
#' @export
normalize.bmop<-function(object,...){
  if (is.null(object$extra$C)){ 
    object$extra$C<-integration_constants(bmop=object)}
  object$ctrpoints<-object$ctrpoint/sum(object$ctrpoints*object$extra$C) 
  return(object)
}



#' Put evidence on a conditional bmop
#' 
#' Normalize a bmop object, makes it integrate to one.
#' @param object a bmop object
#' @param evidence the value of evidence
#' @param evd.pos the position of evidence
#' @param MIN the MIN value as in \code{bmopPar}
#' @param normalize logical, if \code{TRUE} the final bmop object will 
#' be normalized (usually it is not needed since this function is applied to
#' conditional densities)
#' @return a bmop object, the result of imposing some evidence
#' @export
put_evidence.bmop<-function(object,evidence,evd.pos=NULL,
                            MIN=0,normalize=FALSE,...){
 if (length(object$order)==1){return(object)}
 idx<-slice.index(object$ctrpoints,MARGIN = 1)
 if (is.null(evd.pos)){ evd.pos<-2: (length(evidence)+1)}
 bmop<-new_bmop(knots = object$knots[-evd.pos],
                order=object$order[-evd.pos],nk = T)
 bmop$ctrpoints<-apply(object$ctrpoints,MARGIN = -evd.pos,FUN = function(ctr){
   m<-new_bmop(knots=object$knots[evd.pos],order=object$order[evd.pos],nk=T,
               ctrpoints = ctr)
   return(evaluate.bmop(x = evidence,object = m,MIN = MIN))
 })
 if (normalize) {bmop<-normalize.bmop(object = bmop)}
 return(bmop)
}



