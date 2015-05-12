matrixApply <- function(fun,x,y){
  internal <- function(x,y,fun){
    return(apply(y,2,fun,x))
  }
  return(apply(x,2,internal,y,fun))  
}