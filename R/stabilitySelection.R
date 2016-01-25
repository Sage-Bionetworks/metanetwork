stabilitySelection <- function(y,x,alpha=.5,pi=.8,nsamp=100){
  library(glmnet)
  n <- length(y)
  #sampleMatrix <- sample(1:(n*nsamp))
  fxn1 <- function(i,n){
    return(sample(1:n,n))
  }
  sampleMatrix <- sapply(1:n,fxn1,n)
  
  apply(sampleMatrix,2,which)
  
  
  #sampleMatrix <- matrix(sampleMatrix,n,nsamp)
  
  #sampleMatrix <- sampleMatrix%%nsamp+1
  
}