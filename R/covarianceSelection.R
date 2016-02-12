covarianceSelection = function(S,zero){
  library(glasso)
  return(glasso(s=S,rho=0,zero=zero))
}