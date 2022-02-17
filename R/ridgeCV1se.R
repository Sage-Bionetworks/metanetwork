#' Grab CV1se Solution for Ridge Regression
#' 
#' This function runs a glmnet() with ridge regression and pulls the best the beta value of
#' the glmnet.fit object with the largest value of lambda such that error is 
#' within 1 standard error of the minimum.
#' 
#' @inheritParams lassoCV1se 
#' 
#' @return The beta value of the glmnet.fit object with the largest value of 
#' lambda such that error is within 1 standard error of the minimum.
#' @export
ridgeCV1se <- function(y,x,folds=10){
  #require(glmnet)
  res <- glmnet::cv.glmnet(y=y,x=x,nfolds = folds,alpha=0)
  return(res$glmnet.fit$beta[,which(res$lambda==res$lambda.1se)])
}

