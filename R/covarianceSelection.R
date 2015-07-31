covarianceSelection = function(S,zero){
  require(glasso)
  return(glasso(s=S,rho=0,zero=zero))
}