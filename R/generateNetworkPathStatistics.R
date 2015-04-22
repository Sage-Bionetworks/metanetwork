generateNetworkPathStatistics <- function(network,fxn){
  require(dplyr)
  network <- network %>% abs
  cutpoints <- network[network %>% upper.tri %>% which] %>% unique
  internalFunction <- function(cut,network,fxn){
    argList$net <- network>=cut
    (fxn %>% do.call(args=argList)) %>% return
  }
  (cutpoints %>% sapply(internalFunction,network,fxn)) %>% return
}