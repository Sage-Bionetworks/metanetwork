makeModuleNetwork <- function(modules,net){
  #fxn to convert a network into a collapsed module network
  lm <- length(modules)
  moduleNetwork <- matrix(0,lm,lm)
  rownames(moduleNetwork) <- paste0('m',1:lm)
  colnames(moduleNetwork) <- paste0('m',1:lm)
  
  #for (i in 1:lm){
  #  net[modules[[i]],modules[[i]]] <- TRUE
  #}
  for (i in 1:(lm-1)){
    for (j in (i+1):lm){
      if(sum(net[modules[[i]],modules[[j]]])>0){
        moduleNetwork[i,j] <- TRUE
      }
    }
  }
  #diag(net) <- FALSE
  return(moduleNetwork)
}