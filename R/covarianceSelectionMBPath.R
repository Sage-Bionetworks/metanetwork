covarianceSelectionMBPath = function(X,rankedEdges,numberObservations,startI=1){
  #rankedEdges: list of ranked edges
  nedges <- nrow(rankedEdges)
  count <- 1
  bic <- rep(0,nedges)
  cat('bic\n')
  neighborhoods <- vector('list',ncol(X))
  names(neighborhoods)<- colnames(X)
  cat('fastlm\n')
  bicNeighborhood <- apply(X,2,fastlmbic,correction=ncol(X))
  names(bicNeighborhood) <- colnames(X)
  bicCurrent <- sum(bicNeighborhood)
  p <- ncol(X)
  flag <- c()
  while(count<=nedges){
    if(count%%100==0){
      cat('count:',count,'bic:',bicCurrent,'\n')
    }
    gene1 <- colnames(X)[rankedEdges[count,1]]
    gene2 <- colnames(X)[rankedEdges[count,2]]
    neighborhoods[[gene1]] <- c(neighborhoods[[gene1]],gene2)
    neighborhoods[[gene2]] <- c(neighborhoods[[gene2]],gene1)
    
    if(count==startI){
      for(i in 1:ncol(X)){
        bicNeighborhood[i] <- fastlmbic(X[,i],X[,neighborhoods[[i]]],correction=ncol(X))
        #print(i)
      }
      bicCurrent <- sum(bicNeighborhood)
    }
    if(count>startI){
      bicCurrent <- bicCurrent - sum(bicNeighborhood[c(gene1,gene2)])
      bicGene1 <- NA
      bicGene2 <- NA
      try(bicGene1 <- fastlmbic(X[,gene1],X[,neighborhoods[[gene1]]],correction=ncol(X)),silent=T)
      try(bicGene2 <- fastlmbic(X[,gene2],X[,neighborhoods[[gene2]]],correction=ncol(X)),silent=T)
      if(!is.na(bicGene1)){
        bicNeighborhood[gene1] <- bicGene1
      }else{
        flag <- rbind(flag,c(gene1,count))
      }
      if(!is.na(bicGene2)){
        bicNeighborhood[gene2] <- bicGene2
      }else{
        flag <- rbind(flag,c(gene2,count))
      }
      bicCurrent <- bicCurrent + sum(bicNeighborhood[c(gene1,gene2)])
    }
    bic[count] <- bicCurrent
    #plot(bic)
    count <- count+1
  }
  return(list(bic=bic,bicNeighborhood=bicNeighborhood,neighborhoods=neighborhoods,flag=flag))
}