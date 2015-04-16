fastlm <- function(y,X){
  X <- as.matrix(X)
  n1 <- nrow(X)
  X <- cbind(rep(1,n1),X);
  ginv <- solve(t(X)%*%X);
  Xhat <- ginv%*%t(X);
  betahat <- Xhat%*%y;
  sig <- mean((y-X%*%betahat)^2)*((n1)/(n1-ncol(X)));
  zval <- betahat/sqrt(sig*diag(ginv));
  #print('In cleaning')
  return(pt(abs(zval[-1]),n1-2,lower.tail=F)*2);
}