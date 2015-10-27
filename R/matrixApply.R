matrixApply <- function(FUN,MARGIN,X,Y){
  require(dplyr)
  internal <- function(X,Y,FUN,MARGIN){
    return(Y%>% apply(MARGIN,FUN,X))
  }
  return(X %>% apply(MARGIN,internal,Y,FUN,MARGIN))  
}