#' Grab CV1se Solution for Lasso Regression
#' 
#' This function runs a glmnet() with lasso regression and pulls the best the beta value of
#' the glmnet.fit object with the largest value of lambda such that error is 
#' within 1 standard error of the minimum.
#' 
#' @param folds Optional. The number of cross validation folds to partition the 
#' data into. (Default = 10)
#' @inheritParams lassoAIC 
#' 
#' @return The beta value of the glmnet.fit object with the largest value of 
#' lambda such that error is within 1 standard error of the minimum.
#' @export
lassoCV1se <- function(y,x,folds=10){
  #require(glmnet)
  res <- glmnet::cv.glmnet(y=y,x=x,nfolds = folds)
  return(res$glmnet.fit$beta[,which(res$lambda==res$lambda.1se)])
}

#bic:nlogsig + klogn