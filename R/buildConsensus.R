
buildConsensus = function(outputpath,networkFolderId, fileName, pattern_id,bar){
  
  #get all networks from Synapse
  #bar <- synGetFiles(networkFolderId, downloadLocation = outputpath, pattern_id = pattern_id)
  
  loadNetwork <- function(file){
    sparrowNetwork <- data.table::fread(file,stringsAsFactors=FALSE,data.table=F)
    rownames(sparrowNetwork) <- sparrowNetwork$V1
    sparrowNetwork <- sparrowNetwork[,-1]
    gc()
    return(sparrowNetwork)
  }
  getPaths <- function(bar){
    names <- c()
    for (i in 1:length(bar)){
      temp <- bar[[i]]$path
      names <- append(names,temp)
    }
    return(names)
  }
  
  networkFiles <- getPaths(bar)
  networks <- lapply(networkFiles,loadNetwork)
  networks <- lapply(networks,data.matrix)
  
  networks$rankConsensus <- metanetwork::rankConsensus(networks)
  cat('built rank consensus\n')
  cat('write rank consensus\n')
  write.csv(networks$rankConsensus,file=paste0(outputpath,'rankConsensusNetwork.csv'),quote=F)
  if(!is.null(fileName)){
    library(Matrix)
    getNetmethod <- function(networkname){
      temp_names <- strsplit(strsplit(networkname,'Network.csv')[[1]][1],'_')[[1]]
      temp_names <- unlist(temp_names[length(temp_names)])
      return(temp_names)
    }
    networkMethods <- sapply(networkFiles,getNetmethod)
    cat('grabbed methods\n')
    #build rank consensus
    cat('updated methods\n')
    networkMethods <- c(networkMethods,'rankConsensus')
    cat('reading in data\n')
    options(stringsAsFactors = F)
    dataSet <- reader(fileName, row.names=1)
    cat('turning data into data matrix\n')
    dataSet <- data.matrix(dataSet)
    dataSet <- t(dataSet)
    cat('build bicNetworks\n')
    #bicNetworks <- lapply(networks,metanetwork::computeBICcurve,dataSet,maxEdges=1e5)
    bicNetworks <- metanetwork::computeBICcurve(networks$rankConsensus,dataSet,maxEdges=2e5)
    #cat('make names of bicNetworks\n')
    #names(bicNetworks) <- 'rankConsensus'
    cat('save bicNetworks\n')
    save(bicNetworks,file=paste0(outputpath,'bicNetworks.rda'))
  }
  
}
