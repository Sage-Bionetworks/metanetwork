makeModuleNetwork <- function(modules,net){
  #fxn to convert a network into a collapsed module network
  lm <- length(modules)
  for (i in 1:lm){
    net[modules[[i]],modules[[i]]] <- TRUE
  }
  for (i in 1:(lm-1)){
    for (j in i:lm){
      if(sum(net[modules[[i]],modules[[j]]])>0){
        net[modules[[i]],modules[[j]]] <- TRUE
      }
    }
  }
  diag(net) <- FALSE
  return(net)
}