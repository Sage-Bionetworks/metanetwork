###Function to run aracne on the data
aracne <- function(data,path=NULL,pval=NULL,outputpath,na_fill = NULL){
  #path is the a string of the path to th aracne compiled executable
  #data is a matrix of the gene expression data of interest
  installAracne()
  library(caroline)
  if(is.null(pval)){
    pval <- 0.05/choose(nrow(data),2)
  }
  if(is.null(na_fill)){
    data <- na.omit(data)
  }else{
    data[is.na(data)] <- na_fill # Ideally, user could insert a large negative number or use min(data)
  }
  temp_path = paste0(config$input_profile$temp_storage_loc,"/ARACNE")
  setwd(temp_path)

  
  dataMatrix <- cbind(rownames(data),rownames(data),data)
  colnames(dataMatrix) <- c('name1','name2',colnames(data))
  caroline::write.delim(dataMatrix,file='dataMatrix.tsv',sep='\t',quote=FALSE,row.names=FALSE)
  
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
  resultFormat <- lapply(result[-c(1:17)],fun,ref=rownames(data))
  for (i in 1:length(resultFormat)){
    network[resultFormat[[i]]$gene,]<-resultFormat[[i]]$vec;
  }
  #return(network)
  if(pval==1){
    fileName <- paste0(outputpath,'aracneNetwork.csv')
  }else{
    fileName <- paste0(outputpath,'aracneThresholdNetwork.csv')
  }
  #save(network,file=fileName)
  network <- network*upper.tri(network)
  write.csv(network,file=fileName,quote=F)
}


