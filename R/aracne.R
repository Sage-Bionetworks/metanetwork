#' This function applies ARACNE on the data
#' 
#' This function takes in a gene expression matrix of gene expression and applies
#' the ARACNE gene co-expression network analysis frame work to find coexpressed
#' gene pairs in the matrix. The ARACNE framework is also installed from within
#' the package located in the `inst/tools/` direcctory. For more information on the 
#' ARACNE framework see: <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-S1-S7>
#' 
#' @param data Required. The gene expression matrix with rows as sample IDs and
#' columns as Gene or feature IDs 
#' @param path Optional. String containing the path to the aracne compiled 
#' executable. (Default = NULL)
#' @param pval Optional. Cutoff p-value to determine a coexpressed edge. (Default = 0.05)
#' @param outputpath Required. The path the resulting network should be saved to.
#' @param na_fill Optional. Value to replace `NA` values with ideally, a large 
#' negative number or use min(data). (Default = NULL)
#' 
#' @return A Co-Expression Network saved to the path `outputpath` and titled
#'  `aracneThresholdNetwork.csv` (if `pval` < 1) or as `aracneNetwork.csv` if 
#'  `pval` is set to 1. 
#' @export
aracne <- function(data,path=NULL,pval=NULL,outputpath,na_fill = NULL){
  #path is the a string of the path to th aracne compiled executable
  #data is a matrix of the gene expression data of interest
  installAracne()
  #library(caroline)
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


