ridgeBIC <- function(y,x,eigen){
  require(glmnet)
  res <- glmnet(y=y,x=x,alpha=0)
  resid <- y - x%*%res$beta;
  error <- apply(resid^2,2,mean)
  #nonzero <- apply(res$beta,2,function(x) sum(x!=0))
  getDF <- function(lambda,d){
    return(sum(d/(d+lambda)))
  }
  nonzero <- sapply(res$lambda,getDF,eigen)
  n <- nrow(x)
  bic <- n*log(error)+nonzero*2
  
  return(res$beta[,which.min(bic)])
}


