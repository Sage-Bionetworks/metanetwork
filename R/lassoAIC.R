#function to grab BIC solution for lasso
lassoBIC <- function(y,x){
  require(glmnet)
  res <- glmnet(y=y,x=x)
  resid <- y - x%*%res$beta;
  error <- apply(resid^2,2,mean)
  nonzero <- apply(res$beta,2,function(x) sum(x!=0))
  n <- nrow(x)
  bic <- n*log(error)+nonzero*2
  return(res$beta[,which.min(bic)])
}

#bic:nlogsig + klogn