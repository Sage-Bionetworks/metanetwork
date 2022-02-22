#' Fast Linear Modeling BIC 
#' 
#' This function deploys matrix operations to calculate a model BIC given a vector
#' of model coefficients. 
#' 
#' @param x Optional. A numeric vector or matrix of model coefficients. If not set x 
#' becomes a vector of integers the from 1 to length(y). (?)
#' @param y A numeric vector of model coefficients. If not set x 
#' becomes a vector of integers the from 1 to length(y). (?)
#' @param correction A vector of correction factors
#' 
#' @return BIC estimate 
#' @importFrom magrittr %>%
#' @export
fastlmbic <- function(y,x=NULL,correction=1){
  #require(dplyr)
  if(!is.null(x)){
    X <- x %>% as.matrix
    n1 <- X %>% nrow
    X <- (1 %>% rep(n1)) %>% cbind(X)
  }else{
    n1 <- length(y)
    X <- as.matrix(rep(1,n1))
  }
  # Matrix equations to create Beta (coefficients), using X and Y
  ginv <- t(X)%*%X %>% solve();
  Xhat <- ginv%*%t(X);
  betahat <- Xhat%*%y;
  sig <- (((y-X%*%betahat)^2) %>% mean);
  #print('In cleaning')
  # Calculate and return the BIC
  return(n1*(log(sig)+1+log(2*pi))+(ncol(X)+1)*log(n1*correction));
}
