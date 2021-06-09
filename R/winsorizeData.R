# Reference in  https://github.com/Sage-Bionetworks/AMP-AD_Network_Analysis/blob/22bc35fcdfc266fb650bb9b55e89ff70aa22e04c/makeDataForNetworkWG.R
winsorizeData <- function(x){
  library(dplyr)

    winsorize <- function(x,per=.99){
    up <- quantile(x,per,na.rm=T)
    low <- quantile(x,1-per,na.rm=T)
    x[x>=up] <- up
    x[x<=low] <- low
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