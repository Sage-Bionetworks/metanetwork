outerSapply <- function(FUN,X,Y,...){
  require(dplyr)
  internal <- function(X,Y,FUN,...){
    return(Y%>% sapply(FUN,X,...))
  }
  return(X %>% sapply(internal,Y,FUN,...))  
}