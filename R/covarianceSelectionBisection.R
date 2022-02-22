#' Covariance Selection with Bisection Optimization
#'
#' Selects an optimal covariance matrix through BIC convergence.
#'
#' @param numberObservations Required. Number of observations used to calculate 
#' BIC estimates. 
#' @param lowerBoundEdge Required. Numeric specifying the lower bound number of 
#' parameters (d) in BIC calculation: `BIC = -2 * loglikelihood + d * log(N)`
#' @param upperBoundEdge Required. Numeric specifying the upper bound number of 
#' parameters (d) in BIC calculation: `BIC = -2 * loglikelihood + d * log(N)`
#' 
#' @inheritParams covarianceSelection
#' 
#' @return A list containg
#' \itemize{
#' \item{`w` Estimated inverse covariance matrix.}
#'    \item{`resMiddle`	The converged sparse inverse covariance matrix.}
#'    \item{`bicMiddle`	The converged BIC estimate.}
#'    \item{`middleEdge` The converged estimate of parameters.}
#' }
#'
#' @export 
covarianceSelectionBisection = function(S,rankedEdges,numberObservations,lowerBoundEdge,upperBoundEdge){
  #rankedEdges: list of ranked edges
  #library(glasso)
  notConverged <- TRUE
  runStart <- TRUE
  runEnd <- TRUE
  
  getZeros = function(S,rankedEdges,nedges){
    zeros <- which(upper.tri(S)==TRUE,TRUE)
    count <- 1
    while(count <= nedges){
      rowInd <- which(zeros[,1] == rankedEdges[count,1] & zeros[,2] == rankedEdges[count,2])
      #cat(rowInd,'\n')
      zeros <- zeros[-rowInd,]
      count <- count+1
    }
    return(zeros)
  }
  
  while(notConverged){
    if(runStart){
      zerosStart <- getZeros(S,rankedEdges,lowerBoundEdge)
      resStart <- glasso::glasso(S,rho=0,zero=zerosStart)
      bicStart <- -2*resStart$loglik + lowerBoundEdge*log(numberObservations)
    }
    
    middleEdge <- round(lowerBoundEdge/2 + upperBoundEdge/2)
    if(middleEdge == upperBoundEdge){
      notConverged=FALSE
    }
    if(middleEdge == lowerBoundEdge){
      notConverged=FALSE
    }      
    
    cat('middle edge:',middleEdge,'\n')
    zerosMiddle <- getZeros(S,rankedEdges,middleEdge)
    resMiddle <- glasso::glasso(S,rho=0,zero=zerosMiddle)
    bicMiddle <- -2*resMiddle$loglik + middleEdge*log(numberObservations)
    
    if(runEnd){
      zerosEnd <- getZeros(S,rankedEdges,upperBoundEdge)
      resEnd <- glasso::glasso(S,rho=0,zero=zerosEnd)      
      bicEnd <- -2*resEnd$loglik + upperBoundEdge*log(numberObservations)
    }
  
    foo <- which.min((c(bicStart,bicEnd)-bicMiddle))
    
    if(foo==2){
      lowerBoundEdge <- middleEdge
      resStart <- resMiddle
      bicStart <- bicMiddle
      #middleEdge <- round(lowerBoundEdge/2+upperBoundEdge/2)
      runStart <- FALSE
      runEnd <- FALSE
    }else if (foo==1){
      upperBoundEdge <- middleEdge
      resEnd <- resMiddle
      bicEnd <- bicMiddle
      #middleEdge <- round(lowerBoundEdge/2+upperBoundEdge/2)
      runStart <- FALSE
      runEnd <- FALSE 
    }
    
  }
  return(list(res = resMiddle, bic = bicMiddle, nedges = middleEdge))
}