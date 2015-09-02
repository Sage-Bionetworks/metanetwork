#ensembles to networks(small)
ensembleToNetwork <- function(Y){
  internal <- function(i,Y){
    foo <- regressionEnsemble(Y[,i],Y[,-i])
    nullRow <- rep(0,ncol(foo))
    if(i==1){
      return(rbind(nullRow,foo))
    }else if (i>1 & i<ncol(Y)){
      return(rbind(foo[1:(i-1),],nullRow,foo[i:(ncol(Y)-1),]))
    } else{
      return(rbind(foo,nullRow))
    }
  }
  
  
  
  regs <- lapply(1:ncol(Y),internal,Y)
  
  collapseToMatrix <- function(x,i){
    return(x[,i])
  }
  
  biLSapply <- function(FUN,X,Y,...){
    require(dplyr)
    internal <- function(X,Y,FUN,...){
      return(Y%>% sapply(FUN,X,...))
    }
    return(X %>% lapply(internal,Y,FUN,...))  
  }
  model <- list()
  model$matrixRepresentations <- biLSapply(collapseToMatrix,1:ncol(regs[[1]]),regs)
  model$matrixRepresentations <- lapply(model$matrixRepresentations,symmetrisize)
  names(model$matrixRepresentations) <- colnames(regs[[1]])
  require(dplyr)
  w1 <- model$matrixRepresentations[['sparrow2Z']] %>% upper.tri %>% which
  pval <- (model$matrixRepresentations[['sparrow2Z']][w1]^2) %>% pchisq(df=1,lower.tail=F)
  model$nedges <- sum(pval<(0.05/choose(ncol(Y),2)))
  ###rank sum
  foo <- sapply(model$matrixRepresentations,function(x,w1) return(x[w1]),w1)
  foo <- apply(-abs(foo),2,function(x) rank(x,ties.method = 'min'))
  model$aggregateRank <- (rowSums(foo)) %>% rank(ties.method='min')
  model$finalMat <- matrix(0,ncol(Y),ncol(Y))
  model$finalMat[w1] <- model$aggregateRank
  model$finalMat <- symmetrisize(model$finalMat)
  
  return(model)
}
