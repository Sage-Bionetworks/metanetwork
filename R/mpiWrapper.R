mpiWrapper = function(data,nodes,pathv,regressionFunction,outputpath,eigen=NULL,regulatorIndex=NULL,hosts=NULL){
  #initialize MPI
  library('Rmpi');
  #load sparrow library
  library('metaNet');
  nslaves <- nodes;
  #nslaves/nodes: cluster size
  mpi.spawn.Rslaves(nslaves=nslaves,hosts=hosts);
  #if cluster has fewer than 2 nodes, quit

  if (mpi.comm.size() <2){
    cat('More slave processes required.\n');
    mpi.quit();
  }
  
  #clean up function
  .Last <- function(){
    if (is.loaded("mpi_initialize")){
      if (mpi.comm.size(1) > 0){
        cat("Please use mpi.close.Rslaves() to close slaves.\n")
        mpi.close.Rslaves()
      }
      cat("Please use mpi.quit() to quit R\n")
      .Call("mpi_finalize")
    }
  }
  
  foldslave <- function() {
    # Get a task 
    require("metaNet")
    task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
    task_info <- mpi.get.sourcetag() 
    tag <- task_info[2] 
    # While task is not a "done" message. Note the use of the tag to indicate 
    # the type of message
    while (tag != 2) {
      foldNumber <- task$foldNumber
      if(is.null(regulatorIndex)){
        temp_vbsr <- rep(0,p);
        res <- NA;
        set.seed(foldNumber)
        fxnArgs <- list()
        fxnArgs$y <- data[,foldNumber]
        fxnArgs$x <- data[,-foldNumber]
        if(regressionFunction=='vbsrWrapperZ' | regressionFunction=='vbsrWrapperZ2' | regressionFunction=='vbsrWrapperZ2FDR'){
          fxnArgs$n_orderings<-12
        }
        if(!is.null(eigen)){
          fxnArgs$eigen <- eigen
        }
        
        try(res <- do.call(regressionFunction,fxnArgs),silent=TRUE)
        
        if(!is.na(res)){
          temp_vbsr[-foldNumber]<- res;
        }
        temp_res <- c(foldNumber,temp_vbsr);
      }else{
        temp_vbsr <- rep(0,length(regulatorIndex));
        
        if(foldNumber%in%regulatorIndex){
          wi <- which(regulatorIndex%in%foldNumber);
          res <- NA;
          set.seed(foldNumber)
          
          fxnArgs <- list()
          fxnArgs$y <- data[,foldNumber]
          fxnArgs$x <- data[,regulatorIndex][,-wi]
          if(regressionFunction=='sparrowZ' | regressionFunction=='sparrow2Z' | regressionFunction=='sparrow2ZFDR'){
            fxnArgs$n_orderings<-12
          }
          try(res <- do.call(regressionFunction,fxnArgs),silent=TRUE)
          
          if(!is.na(res)){
            temp_vbsr[-wi] <- res;
          }
        }else{
          res <- NA;
          set.seed(foldNumber)
          
          fxnArgs <- list()
          fxnArgs$y <- data[,foldNumber]
          fxnArgs$x <- data[,regulatorIndex]
          if(regressionFunction=='vbsrWrapperZ' | regressionFunction=='vbsrWrapperZ2' | regressionFunction=='vbsrWrapperZ2FDR'){
            fxnArgs$n_orderings<-12
          }
          try(res <- do.call(regressionFunction,fxnArgs),silent=TRUE)
                    
          if(!is.na(res)){
            temp_vbsr <- res;
          }
        }
        temp_res <- c(foldNumber,temp_vbsr);
      }
      
      # Construct and send message back to master
      result <- list(result=temp_res,foldNumber=foldNumber)
      mpi.send.Robj(result,0,1)
      
      # Get a task 
      task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
      task_info <- mpi.get.sourcetag() 
      tag <- task_info[2] 
    }
    
    junk <- 0
    mpi.send.Robj(junk,0,2)
  }
  
  
  p <- ncol(data)
  n <- nrow(data)
  if(regressionFunction=='ridgeAIC'|regressionFunction=='ridgeBIC'){
    eigen <- svd(data)$d^2
  }
  mpi.bcast.Robj2slave(eigen);
  mpi.bcast.Robj2slave(pathv);
  mpi.bcast.Robj2slave(data);
  mpi.bcast.Robj2slave(p);
  mpi.bcast.Robj2slave(n);
  mpi.bcast.Robj2slave(regulatorIndex);
  mpi.bcast.Robj2slave(regressionFunction);

  # Send the function to the slaves
  mpi.bcast.Robj2slave(foldslave)
  # Call the function in all the slaves to get them ready to
  # undertake tasks
  mpi.bcast.cmd(foldslave())
  # Create task list
  tasks <- vector('list')
  for (i in 1:p) {
    tasks[[i]] <- list(foldNumber=i)
  }
  
  # Make the round-robin list for slaves
  n_slaves <- mpi.comm.size()-1
  slave_ids <- rep(1:n_slaves, length=length(tasks))
  
  # Send tasks
  for (i in 1:length(tasks)) {
    slave_id <- slave_ids[i]
    #print(slave_id);
    mpi.send.Robj(tasks[[i]],slave_id,1)
  }
  
  # Collect results
  res_list <- vector("list",p);
  #rssresult <- matrix(0,p,10)
  fold_vec <- 1:p;
  for (i in 1:length(tasks)) {
    message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
    foldNumber <- message$foldNumber
    vv <- which(fold_vec==foldNumber);
    if(length(vv)>0){
      fold_vec <- fold_vec[-vv];
    }
    if((p-length(fold_vec))%%25==0){
      print(p-length(fold_vec));
    }
    results    <- message$result
    res_list[[foldNumber]]<- results
  }
  
  # Perform closing handshake
  for (i in 1:n_slaves) {
    junk <- 0
    mpi.send.Robj(junk,i,2)
  }
  
  for (i in 1:n_slaves) {
    mpi.recv.Robj(mpi.any.source(),2)
  }
  
  # save list to file
  network <- simplify2array(res_list);
  colnames(network) <- colnames(data)
  rownames(network) <- c('fold',colnames(data))
  
  network <- t(network)
  network <- data.frame(network)
  network$fold <- as.integer(network$fold)
  network <- network[,-1]
  if(regressionFunction=='lassoCV1se'|regressionFunction=='lassoCVmin'|regressionFunction=='lassoAIC'|regressionFunction=='lassoBIC'|regressionFunction=='tigress'){
    cat(paste(regressionFunction,sum((network+t(network))!=0)/2,sep=','),'\n',file=paste0(outputpath,'sparsity.csv'),sep='',append=TRUE)
  }else if (regressionFunction=='sparrowZ'){
    #bonferroni
    network2 <- applySparrowBonferroni(network)
    cat(paste('sparrow1Bonferroni',sum(network2!=0)/2,sep=','),'\n',file=paste0(outputpath,'sparsity.csv'),sep='',append=TRUE)
    rm(network2)
    gc()
    network2 <- applySparrowFDR(network)
    cat(paste('sparrow1FDR',sum(network2!=0)/2,sep=','),'\n',file=paste0(outputpath,'sparsity.csv'),sep='',append=TRUE)
    rm(network2)
    gc()
  }else if (regressionFunction=='sparrow2Z'){
    network2 <- applySparrowBonferroni(network)
    cat(paste('sparrow2Bonferroni',sum(network2!=0)/2,sep=','),'\n',file=paste0(outputpath,'sparsity.csv'),sep='',append=TRUE)
    rm(network2)
    gc()
    network2 <- applySparrowFDR(network)
    cat(paste('sparrow2FDR',sum(network2!=0)/2,sep=','),'\n',file=paste0(outputpath,'sparsity.csv'),sep='',append=TRUE)
    rm(network2)
    gc()    
  }
  save(network,file=paste(outputpath,'result_',regressionFunction,'.rda',sep=''));
  
  a <- mpi.close.Rslaves()
  while(a==0){
    a <- mpi.close.Rslaves()
  }
  #mpi.bcast.cmd(q("no"));
  mpi.quit(save="no")
  
}