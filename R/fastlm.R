#' Fast Linear Modeling
#' 
#' This function retuns results from a fast linear model
#' @param y Required. A vector of response values 
#' @param x Required. A vector or matric of the same number of observations/rows
#' as y
#' @return A vector of ZValues
#' 
#' @importFrom magrittr %>%
#' @export
fastlm <- function(y,x){
  #require(dplyr)
  X <- x %>% as.matrix
  n1 <- X %>% nrow
  X <- (1 %>% rep(n1)) %>% cbind(X)
  ginv <- t(X)%*%X %>% solve();
  Xhat <- ginv%*%t(X);
  betahat <- Xhat%*%y;
  sig <- (((y-X%*%betahat)^2) %>% mean)*((n1)/(n1- X %>% ncol));
  zval <- betahat/((sig*(ginv %>% diag)) %>% sqrt);
  #print('In cleaning')
  return(zval[-1]);
}