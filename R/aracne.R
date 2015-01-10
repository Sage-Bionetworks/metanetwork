###Function to run aracne on the data
aracne <- function(data,path=NULL,pval=NULL){
  #path is the a string of the path to th aracne compiled executable
  #data is a matrix of the gene expression data of interest
  if(!is.null(path)){
    setwd(path)
  }
  if(is.null(pval)){
    pval <- 0.05/choose(ncol(data),2)
  }
  
  dataMatrix <- cbind(rownames(data),rownames(data),data)
  colnames(dataMatrix) <- c('name1','name2',colnames(data))
  write.table(dataMatrix,file='dataMatrix.tsv',sep='\t',quote=FALSE,row.names=FALSE)
  
  command_string <- paste('./aracne2 -i dataMatrix.tsv -a adaptive_partitioning -p ',pval,' -o result.out',sep='')
  system(command_string) 
  result <- readLines('result.out')
  network <- matrix(0,nrow(data),nrow(data))
  rownames(network) <- rownames(data)
  colnames(network) <- rownames(data)
  fun <- function(x,ref){
    vec <- strsplit(x,'\t')[[1]]
    model <-list();
    model$gene <- vec[1];
    model$vec <- rep(0,length(ref));
    names(model$vec) <- ref;
    model$keepGene <- vec[-1][1:length(vec[-1])%%2==1];
    model$weights <- vec[-1][1:length(vec[-1])%%2==0];
    model$vec[model$keepGene] <- model$weights
    model$vec <- as.numeric(model$vec)
    return(model)
  }
  resultFormat <- lapply(result[-c(1:18)],fun,ref=rownames(data))
  for (i in 1:length(resultFormat)){
    network[resultFormat[[i]]$gene,]<-resultFormat[[i]]$vec;
  }
  return(network)
}

n <- 20
m <- 100

dum <- matrix(rnorm(n*m),n,m)
colnames(dum) <- paste('a',1:m,sep='')
rownames(dum) <- paste('b',1:n,sep='')
#dum <- cbind(rownames(dum),dum)
a<-aracne(data=dum,pval = 0.1)
write.table(dum,file='dum.txt',sep='\t',quote=FALSE)
system('./aracne2 -i dum.txt -a adaptive_partitioning -p 1e-3')
system('./aracne2 -i test/arraydata10x336.exp -a adaptive_partitioning -p 1e-3')
