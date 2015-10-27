regressionEnsemble <- function(y,
                               x,
                               methods=c('ridgeAIC','ridgeBIC','ridgeCV1se','ridgeCVmin','genie3','lassoCVmin','lassoCV1se','tigressRootN','sparrowZ','sparrow2Z','lassoAIC','lassoBIC','elasticNetCVmin','elasticNetCV1se','elasticNetAIC','elasticNetBIC'),eigen=NULL){
  internal <- function(fxn,y,x,eigen=NULL){
    if(fxn=='ridgeAIC'|fxn=='ridgeBIC'){
      eigen <- svd(x)$d^2
    }
    if(is.null(eigen)){
      return(do.call(fxn,list(y=y,x=x)))
    }else{
      return(do.call(fxn,list(y=y,x=x,eigen=eigen)))
    }
  }
  modelParams <- (sapply(methods,internal,y,x))
  colnames(modelParams)
  return(modelParams)
}