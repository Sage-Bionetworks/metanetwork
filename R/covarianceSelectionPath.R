covarianceSelectionPath = function(S,rankedEdges,numberObservations){
  #rankedEdges: list of ranked edges
  library(glasso)
  zeros <- which(upper.tri(S)==TRUE,TRUE)
  count <- 1
  bic <- nrow(rankedEdges)
  
  
  while(count<=nrow(rankedEdges)){
    rowInd <- which(zeros[,1]==rankedEdges[count,1] & zeros[,2] == rankedEdges[count,2])
    cat(rowInd,'\n')
    zeros <- zeros[-rowInd,]
    if(count==1){
      res <- glasso(S,rho=0,zero=zeros)
    }else{
      res <- glasso(S,rho=0,zero=zeros,w.init = res$w, wi.init = res$wi,start = 'warm')
    }
    bic[count] <- -2*res$loglik + count*log(numberObservations)
    count <- count+1
  }
  
  
  bicmin <- min(bic)
  nedges <- which.min(bic)
  zeros <- which(upper.tri(S)==TRUE,TRUE)
  count <- 1
  while(count<=nedges){
    rowInd <- which(zeros[,1]==rankedEdges[count,1] & zeros[,2] == rankedEdges[count,2])
    cat(rowInd,'\n')
    zeros <- zeros[-rowInd,]
    count <- count+1
  }
  res <- glasso(S,rho=0,zero=zeros,w.init = res$w, wi.init = res$wi,start = 'warm')
  
  return(list(bic=bic,bicmin=bicmin,w=res$w,wi=res$wi))
}