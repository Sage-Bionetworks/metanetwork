#' Grab CV1min Solution for Lasso Regression
#' 
#' This function runs a glmnet() with lasso regression and pulls the best the beta 
#' value of the glmnet.fit object with value of lambda that gives minimum cvm.
#' 
#' @inheritParams lassoCV1se 
#' 
#' @return The beta value of the glmnet.fit object with value of lambda that gives 
#' minimum cvm
#' @export
lassoCVmin <- function(y,x,folds=10){
  #require(glmnet)
  res <- glmnet::cv.glmnet(y=y,x=x,nfolds = folds)
  return(res$glmnet.fit$beta[,which(res$lambda==res$lambda.min)])
}
