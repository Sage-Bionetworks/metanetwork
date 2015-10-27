metaReg <- function(y,x,eigen){
  print(system.time(lassoRes <- lassoBIC(y=y,x=x)))
  print(system.time(ridgeRes <- ridgeBIC(y=y,x=x,eigen=eigen)))
  print(system.time(vbsrRes <- vbsrWrapper(y=y,x=x)))
  print(system.time(rfRes <- genie3(y=y,x=x)))
  print(system.time(ssRes <- tigress(y=y,x=x)))
  return(cbind(vbsrRes,lassoRes,ssRes,ridgeRes,rfRes))
}
