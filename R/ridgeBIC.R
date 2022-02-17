#' Grab BIC Solution for Ridge Regression
#' 
#' This function runs a glmnet() function with ridge regression and pulls the best 
#' BIC estimate.
#' 
#' @inheritParams lassoAIC 
#' 
#' @return Lowest BIC value.
#' 
#' @export
ridgeBIC <- function(y,x,eigen){
  #require(glmnet)
  res <- glmnet::glmnet(y=y,x=x,alpha=0)
  resid <- y - x%*%res$beta;
  error <- apply(resid^2,2,mean)
  #nonzero <- apply(res$beta,2,function(x) sum(x!=0))
  getDF <- function(lambda,d){
    return(sum(d/(d+lambda)))
  }
  nonzero <- sapply(res$lambda,getDF,eigen)
  n <- nrow(x)
  bic <- n*log(error)+nonzero*log(n)
  
  return(res$beta[,which.min(bic)])
}


