
###function to run sparrow on the data
sparrow <- function(data,driverIndex,...){
  extractEdges <- function(index,data,driverIndex,...){
    X <- data;
    if(index%%100==0){
      cat(index,'!','\n')
    }
    if (index%in%driverIndex){
      y <- X[,index];
      G <- X[,driverIndex];
      wi <- which(driverIndex%in%index);
      G <- G[,-wi];
      res <- vbsr(y=y,X=G,...);
      pval <- rep(1,length(res$pval)+1);
      pval[-wi] <- res$pval;
      beta <- rep(0,length(res$beta)+1);
      beta[-wi] <- res$beta;
      
    }else{
      y <- X[,index];
      G <- X[,driverIndex];
      res <- vbsr(y=y,X=G,...);
      pval <- res$pval;
      beta <- res$beta;
    }
    return(list(pval=pval,beta=beta));
  }
  ind <- 1:ncol(data);
  edgeMat <- lapply(ind,extractEdges,data=data,driverIndex=driverIndex,...)
  #print(edgeMat)
  edgeMatPval <- sapply(edgeMat,function(x){return(x$pval)});
  edgeMatBeta <- sapply(edgeMat,function(x){return(x$beta)});
  driverMat <- edgeMatPval< 0.05/(ncol(data)*length(driverIndex));
  colnames(driverMat) <- colnames(data);
  rownames(driverMat) <- names(driverIndex)
  colnames(edgeMatBeta) <- colnames(driverMat)
  rownames(edgeMatBeta) <- rownames(driverMat)
  return(list(driverMat=driverMat,edgeMatBeta=edgeMatBeta))
}