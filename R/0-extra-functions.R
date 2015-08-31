
removeNA<-function(data){
  return(as.data.frame(data[!apply(X = data,MARGIN = 1,function(d){
    any(is.na(d))
  }),]))
}




multiplicity<-function(x,v){
  return(length(v[x==c(x,v)])-1)
}


locate<-function(a,v){
  v<-sort(c(a,v))
  return(max((1:(length(v)))[v==a]-1))
}



fix_data<-function(data){
 return(as.data.frame(data))
}


#' Bins grouping
#' 
#' This function extend the histogram class for multi-dimensional datasets
#' @param data a dataset
#' @param breaks function or positive integer
#' @param ... additional parameters
#' @return an object \code{bins}, the data are grouped into bins uniformly
#' @export
as.bins<-function(data,breaks=nclass.FD,...){
  bins<-list()
  data<-as.data.frame(data)
  if (dim(data)[2]==1){
    data<-as.matrix(data)
    dim(data)<-NULL
    H<-hist(data,breaks=breaks,plot = F)
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
  Ns<-sapply(Mids,length)
  bins$mids<-expand.grid(Mids)
  bins$counts<-c(rep(0,dim(bins$mids)[1]))
  Es<-c()
  Es[1]<-1
  for (i in 2:(length(Ns))){
    Es[i]<-Es[i-1]*Ns[i-1]
  }
  for (i in 1:(dim(data)[1])){
    pos<-sapply(1:(dim(data)[2]),FUN = function(j){
      locate(data[i,j],Seqs[[j]])
    })
    pos<-pos-1
    pos[1]<-pos[1]+1
    idx<-sum( (pos)*Es )
    bins$counts[idx] <- bins$counts[idx] + 1
  }
  idx<-bins$counts!=0
  bins$counts<-bins$counts[idx]
  bins$mids<-bins$mids[idx,]
  class(bins)<-"bins"
  return(bins)
  
}



