covarianceSelectionPath = function(S,rankedEdges){
  #rankedEdges: list of ranked edges
  library(glasso)
  zeros <- which(upper.tri(S)==TRUE,TRUE)
  count <- 1
  bic <- nrow(rankedEdges)
  while(count<=nrow(rankedEdges)){
    rowInd <- which(zeros[,1]==rankedEdges[count,1] & zeros[,2] == rankedEdges[count,2])
    zeros <- zeros[-rowInd,]
    if(count==1){
      res <- glasso(S,rho=0,zeros=zeros)
    }else{
      
    }
    
    
    count <- count+1
  }
  
}