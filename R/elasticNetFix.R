elasticNetFix <- function(y,x,lambdaGrid=exp(seq(-4,4,length.out=100)),alphaGrid=seq(0,1,length.out=10)){
  #foldid <- sample(1:10,length(y),replace = T)
  foldid <- c(sapply(1:10,function(x,y){return(rep(x,y))},ceiling(length(y)/10)))
  foldid <- foldid[sample(1:length(foldid))]
  foldid <- foldid[1:length(y)]
  library(glmnet)
  elasticNetWrapFxn = function(alpha,y,x,lambda,foldid){
    res <- cv.glmnet(y=y,x=x,alpha=alpha,lambda=lambda,foldid=foldid)
    k1 <- which(res$lambda==res$lambda.1se)
    return(res$cvm)
  }
  foo <- sapply(alphaGrid,elasticNetWrapFxn,y,x,lambdaGrid,foldid)
  foo2 <- which.min(foo)
  res <- cv.glmnet(y=y,x=x,alpha=alphaGrid[foo2])
  return(res)
}