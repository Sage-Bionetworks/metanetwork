stabilitySelection <- function(y,x,alpha=.5,pi=.8,nsamp=100){
  library(glmnet)
  n <- length(y)
  sampleMatrix <- sample(1:(n*nsamp))
  #sampleMatrix <- matrix(sampleMatrix,n,nsamp)
  
  #sampleMatrix <- sampleMatrix%%nsamp+1
  
}