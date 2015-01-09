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
  
  
  
}

n <- 100
m <- 100

dum <- matrix(rnorm(n*m),n,m)
colnames(dum) <- paste('a',1:100,sep='')
rownames(dum) <- paste('b',1:100,sep='')
#dum <- cbind(rownames(dum),dum)
aracne(data=dum)
write.table(dum,file='dum.txt',sep='\t',quote=FALSE)
system('./aracne2 -i dum.txt -a adaptive_partitioning -p 1e-3')
system('./aracne2 -i test/arraydata10x336.exp -a adaptive_partitioning -p 1e-3')
