#' Grab BIC Solution for Lasso
#' 
#' This function runs a glmnet() function with lasso regression and pulls the best 
#' BIC estimate.
#' 
#' @inheritParams lassoAIC 
#' 
#' @return Lowest BIC value.
#' 
#' @export
lassoBIC <- function(y,x){
  #require(glmnet)
  res <- glmnet::glmnet(y=y,x=x)
  resid <- y - x%*%res$beta;
  error <- apply(resid^2,2,mean)
  nonzero <- apply(res$beta,2,function(x) sum(x!=0))
  n <- nrow(x)
  bic <- n*log(error)+nonzero*log(n)
  return(res$beta[,which.min(bic)])
}