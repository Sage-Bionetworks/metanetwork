#function to grab BIC solution for lasso
lassoBIC <- function(y,X){
  require(glmnet)
  res <- glmnet(y=y,x=X)
  resid <- y - X%*%res$beta;
  error <- apply(resid^2,2,mean)
  nonzero <- apply(res$beta,2,function(x) sum(x!=0))
  n <- nrow(X)
  bic <- n*log(error)+nonzero*log(n)
  return(res$beta[,which.min(bic)])
}

#bic:nlogsig + klogn