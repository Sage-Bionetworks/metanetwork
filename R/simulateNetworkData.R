#' Simulate Network Data
#' 
#' This function simulates a gene correlation network and inverse covariance matrix
#' of user determined size.
#' 
#' @param n Required. The number of rows to simulate in the data. These correspond 
#' to the number of simulated samples
#' @param p Required. The number of columns to simulate in the data. These correspond 
#' to the number of gene features
#' @param prop Required. a value between 0 and 1, corresponding the probability
#' value in the `rbinom()` adjustment of expression values.
#' @param adjustment Optional. A matrix, vector, or 1D array, or missing. This
#' corresponds to the x vaule in the `diag()` adjustment of the inverse covariance
#' matrix. (Default = 10)
#' 
#' @return A list containing a simulated network and the inverse covariance matrix
#' of the simulated network.s
#' 
#' @export
simulateNetworkData <- function(n,p,prop,adjustment=1e1){
  #library(MASS)
  
  LAM <- matrix(stats::rnorm(p^2)*stats::rbinom(p^2,1,prop),p,p)
  diag(LAM) <- 1;
  OM <- LAM%*%t(LAM)
  OM <- OM + diag(adjustment,p)
  
  
  
  Y <- MASS::mvrnorm(n,rep(0,p),solve(OM))
  
  rn <- sapply(1:n,function(i){return(paste(sample(letters,5,replace=T),collapse=''))})
  cn <- sapply(1:p,function(i){return(paste(sample(letters,5,replace=T),collapse=''))})
  
  while(sum(duplicated(rn))>0){
    k <- sum(duplicated(rn))
    rn[duplicated(rn)] <- sapply(1:k,function(i){return(paste(sample(letters,5,replace=T),collapse=''))})
  }
  
  while(sum(duplicated(cn))>0){
    k <- sum(duplicated(cn))
    cn[duplicated(cn)] <- sapply(1:k,function(i){return(paste(sample(letters,5,replace=T),collapse=''))})
  }
  
  rownames(Y) <- rn
  colnames(Y) <- cn
  
  rownames(OM) <- cn
  colnames(OM) <- cn
  
  return(list(data=Y,inverseCovariance=OM))
  
}