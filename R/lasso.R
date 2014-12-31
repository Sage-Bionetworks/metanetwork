###Function to run lasso on the data
lasso <- function(data,driverIndex,...){
  require(glmnet)
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
      res <- lassoBIC(y=y,X=G,...);
      #pval <- rep(1,length(res$pval)+1);
      #pval[-wi] <- res$pval;
      beta <- rep(0,length(res)+1);
      beta[-wi] <- res;
      
    }else{
      y <- X[,index];
      G <- X[,driverIndex];
      res <- vbsr(y=y,X=G,...);
      #pval <- res$pval;
      beta <- res;
    }
    return(list(beta=beta));
  }
  ind <- 1:ncol(data);
  edgeMat <- lapply(ind,extractEdges,data=data,driverIndex=driverIndex,...);
  edgeMatBeta <- sapply(edgeMat,function(x){return(x$beta)});
  return(edgeMatBeta);
}