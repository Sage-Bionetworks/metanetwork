c3netWrapper  = function(data){
  library(c3net)
  network <- c3net(t(data))
  save(network,file='result_c3net.rda')
}