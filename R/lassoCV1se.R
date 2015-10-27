#function to grab BIC solution for lasso
lassoCV1se <- function(y,x,folds=10){
  require(glmnet)
  res <- cv.glmnet(y=y,x=x,nfolds = folds)
  return(res$glmnet.fit$beta[,which(res$lambda==res$lambda.1se)])
}

#bic:nlogsig + klogn