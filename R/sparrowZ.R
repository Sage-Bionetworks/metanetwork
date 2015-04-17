#vbsrWrapper
sparrowZ <- function(y,x,...){
  require(vbsr)
  result <- vbsr(y=y,X=x,...)$z
  return(result)
}
