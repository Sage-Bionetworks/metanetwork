fastlm <- function(y,x){
  require(dplyr)
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