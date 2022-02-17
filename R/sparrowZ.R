#' Runs vbsr gene across a matrix 
#' 
#' This function wraps variable bays spike regression of a genes expression across
#' a matrix of genes expressed in the same samples.
#' 
#' @param y Required. response variable. Normally distributed errors for 
#' family="normal". For family="binomial" should be coded as a vector of 0's and 1's.
#' @param x Required. Design matrix, an n x m matrix, with rows as observations.
#'
#' @return A coexpression value
#' @export
sparrowZ <- function(y,x,...){
  require(vbsr)
  result <- vbsr(y=y,X=x,...)$z
  return(result)
}
