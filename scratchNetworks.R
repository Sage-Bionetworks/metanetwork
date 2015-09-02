#n <- as.numeric(commandArgs(TRUE)[[1]])
#p <- as.numeric(commandArgs(TRUE)[[2]])
#prop <- as.numeric(commandArgs(TRUE)[[3]])
#file <- as.character(commandArgs(TRUE)[[4]])



p <- 25
prop <- 2/25

require(MASS)
LAM <- matrix(rnorm(p^2)*rbinom(p^2,1,prop),p,p)
diag(LAM) <- 1;
OM <- LAM%*%t(LAM)
OM <- OM + diag(1,p)

aggregateRank = function(x){
  rankMat <- apply(-x,2,rank,ties.method='min')
  return(rank(rowMeans(rankMat),ties.method='min'))
}

n <- 200
Y <- mvrnorm(n,rep(0,p),solve(OM))
Y <- Y+rnorm(n)
result <- ensembleToNetwork(Y)

res <- regressionEnsemble(Y[,1],Y[,-1])
res <- abs(res)
res <- cbind(res,rowSums(scale(res)))
res <- cbind(apply(-res,2,rank),aggregateRank(res[,-ncol(res)]))
plot(as.factor((OM[-1,1])!=0),res[,7])
plot(as.factor((OM[-1,1])!=0),res[,ncol(res)])
#cor(cbind(scale(res),OM[-1,1]!=0),method = 'spearman')
#pairs(cbind(scale((res)),OM[-1,1]!=0))

library(ROCR)
#getAUC <- function(pred,lab){
#  return(performance(prediction(-pred,lab),'auc')@y.values[[1]])
#}
#auc <- apply(res,2,getAUC,OM[-1,1]!=0)
#(auc-.5)*2
pred1 <- prediction(-abs(result$finalMat[w1]),labels=OM[w1]!=0)
#pred1 <- prediction(abs(result$matrixRepresentations[[16]][w1]),labels=OM[w1]!=0)
perf1 <- performance(pred1,'auc')
perf2 <- performance(pred1,'prec','rec')
cat(perf1@y.values[[1]],'\n')
plot(perf2)
plot(as.factor(OM[w1]!=0),result$aggregateRank)
