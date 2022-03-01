#' Winsorize a Gene Expression Matrix
#' 
#' Will winsorize a given gene expression matrix
#' 
#' @param x Gene expression matrix 
#' 
#' @return A winsorized expression matrix
#' 
#' @export
winsorizeData <- function(x){
  winsorize <- function(x,per=.99){
    up <- stats::quantile(x,per,na.rm=T)
    low <- stats::quantile(x,1-per,na.rm=T)
    x[x >= up] <- up
    x[x <= low] <- low
    return(x)
  }

  replaceNaMean <- function(x){
    if(sum(is.na(x))>0){
      y <- x
      y[is.na(x)] <- mean(x,na.rm=T)
      return(y)
    }else{
      return(x)
    }
  }
  
  x <- t(x)
  x <- apply(x,2,winsorize)
  x <- apply(x,2,replaceNaMean)
  x <- scale(x)
  return(x)
}