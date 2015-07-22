makeModuleNetwork <- function(modules,net){
  #fxn to convert a network into a collapsed module network
  for (i in 1:length(modules)){
    net[modules[[i]],modules[[i]]] <- TRUE
  }
  diag(net) <- FALSE
  return(net)
}