sparrowMPI = function(data,nodes,pathv,regulatorIndex=NULL,hosts=NULL){
  #initialize MPI
  library('Rmpi');
  #load sparrow library
  library('vbsr');
  nslaves <- nodes;
  #nslaves/nodes: cluster size
  mpi.spawn.Rslaves(nslaves=nslaves,hosts=hosts);
  #cat('well at least we made it here!')
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
    require("vbsr")
    task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
    task_info <- mpi.get.sourcetag() 
    tag <- task_info[2] 
    #cat('in fold slave tag:',tag,'\n')
    # While task is not a "done" message. Note the use of the tag to indicate 
    # the type of message
    while (tag != 2) {
      # Perform the task.  
      foldNumber <- task$foldNumber
      #rss <- double(p)
      #for (i in 1:p) {
      #cat('in fold',foldNumber,'\n')
      if(is.null(regulatorIndex)){
        #cat('running vbsr\n')
        temp_vbsr <- rep(0,p);
        
        #temp_cor <- rep(0,p);
        #set.seed(1);
        res <- NA;
        set.seed(foldNumber)
        #print(data[,foldNumber])
        #print(data[,-foldNumber][1:5,1:5])
        #print(dim(data[,-foldNumber]))
        res <- vbsr(y=data[,foldNumber],X=data[,-foldNumber],n_orderings=12);
        #cat('The system works\n')
        if(!is.na(res)){
          temp_vbsr[-foldNumber]<- res$z;
          #temp_cor[-foldNumber] <- res$cor;
        }
        temp_res <- c(foldNumber,temp_vbsr);
        #fn <- paste(pathv,'gene',foldNumber,'.out',sep='')
        #cat(temp_res,file=fn);
      }else{
        #temp_res <- rep(0,length(regulatorIndex));
        temp_vbsr <- rep(0,length(regulatorIndex));
        #temp_cor <- rep(0,length(regulatorIndex));
        
        if(foldNumber%in%regulatorIndex){
          wi <- which(regulatorIndex%in%foldNumber);
          #set.seed(1);
          res <- NA;
          set.seed(foldNumber)
          try(res <- vbsr(data[,foldNumber],data[,regulatorIndex][,-wi],n_orderings=12),silent=TRUE)
          if(!is.na(res)){
            temp_vbsr[-wi] <- res$z;
            #temp_cor[-wi] <- res$cor;
          }
        }else{
          #set.seed(1);
          res <- NA;
          set.seed(foldNumber)
          try(res <- vbsr(data[,foldNumber],data[,regulatorIndex],n_orderings=12),silent=TRUE);
          if(!is.na(res)){
            temp_vbsr <- res$z;
            #temp_cor <- res$cor;
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
  
  
  #load("test_data.RData");
  
  #n <- nrow(Y);
  #p <- ncol(Y);
  #p <- 100;
  
  
  # Now, send the data to the slaves
  p <- ncol(data)
  n <- nrow(data)
  mpi.bcast.Robj2slave(pathv);
  mpi.bcast.Robj2slave(data);
  mpi.bcast.Robj2slave(p);
  mpi.bcast.Robj2slave(n);
  mpi.bcast.Robj2slave(regulatorIndex);
  #cat('data was sent to slaves\n')
  #mpi.bcast.Robj2slave(regulatorIndex);
  #mpi.bcast.Robj2slave(thedata)
  #mpi.bcast.Robj2slave(fold)
  #mpi.bcast.Robj2slave(p)
  
  # Send the function to the slaves
  mpi.bcast.Robj2slave(foldslave)
  #cat('foldslave was sent to slaves\n')
  # Call the function in all the slaves to get them ready to
  # undertake tasks
  mpi.bcast.cmd(foldslave())
  #cat('foldslave was called in all slaves\n')
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
    print(p-length(fold_vec));
    results    <- message$result
    #cat(results[1:5],'\n')
    #rssresult[,foldNumber] <- results
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
  result <- simplify2array(res_list);
  #result <- result[,order(result[1,])]
  colnames(result) <- paste0(colnames(data),'_dep')
  rownames(result) <- c('fold',paste0(colnames(data),'_indep'))
  result <- t(result)
  result <- data.frame(result)
  result$fold <- as.integer(result$fold)
  save(result,file=paste(pathv,'result.rda',sep=''));
  
  mpi.close.Rslaves()
  #mpi.bcast.cmd(q("no"));
  mpi.quit(save="no")
  
}