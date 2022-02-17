#' Covariance Selection with Graphical Lasso
#'
#' Estimates a sparse inverse covariance matrix using a lasso (L1) penalty.
#'
#' @param S Required. A symetric p-by-p covariance matrix.
#' @param rankedEdges Required. A list of ranked edges to be constrained by zero.
#' 
#' @return A list with components.
#' \itemize{
#'    \item{`w` Estimated inverse covariance matrix.}
#'    \item{`loglik`	Value of maximized log-likelihodo+penalty.}
#'    \item{`errflag`	Memory allocation error flag: 0 means no error; !=0 means 
#'    memory allocation error - no output returned.}
#'    \item{`approx`	Value of input argument approx.}
#'    \item{`del` Change in parameter value at convergence.}
#'    \item{`niter`	Number of iterations of outer loop used by algorithm.}
#' }
#'
#' @export 
covarianceSelection = function(S,rankedEdges){
  #library(glasso)
  
  nedges <- nrow(rankedEdges)
  zeros <- which(upper.tri(S)==TRUE,TRUE)
  count <- 1
  while(count<=nedges){
    rowInd <- which(zeros[,1]==rankedEdges[count,1] & zeros[,2] == rankedEdges[count,2])
    #cat(rowInd,'\n')
    zeros <- zeros[-rowInd,]
    count <- count+1
  }
  res <- glasso::glasso(S,rho=0,zero=zeros)
  return(res)
}