fastlm <- function(y,x){
  require(dplyr)
  X <- as.matrix(x)
  n1 <- nrow(X)
  X <- cbind(rep(1,n1),X);
  ginv <- t(X)%*%X %>% solve();
  Xhat <- ginv%*%t(X);
  betahat <- Xhat%*%y;
  sig <- mean((y-X%*%betahat)^2)*((n1)/(n1-ncol(X)));
  zval <- betahat/((sig*diag(ginv)) %>% sqrt());
  #print('In cleaning')
  return(zval[-1]);
}