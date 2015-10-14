#vbsrWrapper
vbsrWrapper <- function(y,x,...){
  require(vbsr)
  result <- vbsr(y=y,X=x,...)$post
  return(result)
}
