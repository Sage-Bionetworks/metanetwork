covarianceSelection = function(S,rankedEdges){
  library(glasso)
  
  nedges <- nrow(rankedEdges)
  zeros <- which(upper.tri(S)==TRUE,TRUE)
  count <- 1
  while(count<=nedges){
    rowInd <- which(zeros[,1]==rankedEdges[count,1] & zeros[,2] == rankedEdges[count,2])
    #cat(rowInd,'\n')
    zeros <- zeros[-rowInd,]
    count <- count+1
  }
  res <- glasso(S,rho=0,zero=zeros)
  return(res)
}