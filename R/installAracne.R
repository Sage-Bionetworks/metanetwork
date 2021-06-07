installAracne <- function(path=NULL){
  arc_path = system.file("tools", "ARACNE.src.tar.gz", package = "metanetwork")
  system(paste('tar -xzvf ',arc_path,sep=''))
  #str1 <- aracne@filePath;
  #str2 <- strsplit(str1,'ARACNE.src.tar.gz')[[1]]
  str3 <- 'ARACNE/'
  cat(str3,'\n')  
  setwd(str3)
  #system(paste(str2,'/ARACNE.src/ARACNE/',sep=''))
  system('make')
  #setwd('../')
  #return(paste(getwd(),'/ARACNE/',sep=''))
}
