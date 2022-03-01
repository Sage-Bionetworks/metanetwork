#' Grab CV1min Solution for Ridge Regression
#' 
#' This function runs a glmnet() with ridge regression and pulls the best the beta 
#' value of the glmnet.fit object with value of lambda that gives minimum cvm.
#' 
#' @inheritParams lassoCVmin 
#' 
#' @return The beta value of the glmnet.fit object with value of lambda that gives 
#' minimum cvm
#' @export
ridgeCVmin <- function(y,x,folds=10){
  #require(glmnet)
  res <- glmnet::cv.glmnet(y=y,x=x,nfolds = folds,alpha=0)
  return(res$glmnet.fit$beta[,which(res$lambda==res$lambda.min)])
}
