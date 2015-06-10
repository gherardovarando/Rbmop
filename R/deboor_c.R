



#############################################################
#                                                           
#   deboor_c function                                  
#   Author: Gherardo Varando, Luis Rodriguez-lujan, Donald-duck
#   Date: April, 30, 2015                                  
#   Version: 0.1                                          
#                                                           
#############################################################
#' @useDynLib Rbmop deboor_eval
deboor_c <- function(t,k,knots,ctr,MIN=1E-10){

    l_knots <- length(knots)
    l_ctr <- length(ctr)
    
    res <- .C("deboor_eval",t = as.double(t),
              k = as.integer(k),
              knots_len = as.integer(l_knots),
              knots = as.double(knots),
              ctr_len = as.integer(l_ctr),
              ctr = as.double(ctr),
              min = as.double(MIN),
              retval = as.double(0))
    return(res$retval)
}

