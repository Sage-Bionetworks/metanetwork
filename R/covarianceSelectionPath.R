#' Covariance Selection with Path BIC Selection
#'
#' Selects an optimal covariance matrix through Path BIC Selection.
#
#' @param S Required. An expression matrix
#' @param rankedEdges Required. An edge list in the form of a matrix where column
#' one is gene one and column 2 is gene two
#' @param numberObservations Required. The number of samples comprising the network (?)
#' @param n Optional. Start at the first edge in `rankedEdges` (Default = 1)
#' 
#' @return A list object of containing the BIC estimate, bicNeighborhood (?) , neighborhoods (?), flag (?)
#' @export A list containg
#' \itemize{
#'     \item{`bic` Estimated BIC.}
#'    \item{`bicmin` Lowest Estimated BIC.}
#'    \item{`w`	The covariance matrix.}
#'    \item{`wi` The inverse covariance matrix.}
#' }
covarianceSelectionPath = function(S,rankedEdges,numberObservations,n){
  #rankedEdges: list of ranked edges
  #library(glasso)
  zeros <- which(upper.tri(S)==TRUE,TRUE)
  count <- 1
  bic <- nrow(rankedEdges)
  
  
  while(count<=nrow(rankedEdges)){
    rowInd <- which(zeros[,1]==rankedEdges[count,1] & zeros[,2] == rankedEdges[count,2])
    #cat(rowInd,'\n')
    zeros <- zeros[-rowInd,]
    if(count==1){
      res <- glasso::glasso(S,rho=0,zero=zeros)
    }else{
      res <- glasso::glasso(S,rho=0,zero=zeros,w.init = res$w, wi.init = res$wi,start = 'warm')
    }
    print(res$w[1:5,1:5])
    print(res$wi[1:5,1:5])
    #bic[count] <- -2*res$loglik + count*log(numberObservations)
    loglik <- (-n/2)*(sum(log(eigen(res$w)$values)) + sum(diag(S%*%res$wi))+nrow(S)*(2*pi))
    bic[count] <- -2*(loglik) + count*log(numberObservations)
    cat(count,bic[count],'\n')
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
  res <- glasso::glasso(S,rho=0,zero=zeros,w.init = res$w, wi.init = res$wi,start = 'warm')
  
  return(list(bic=bic,bicmin=bicmin,w=res$w,wi=res$wi))
}