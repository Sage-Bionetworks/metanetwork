installAracne <- function(path=NULL){
  cwd <- getwd()
  aracne <- synGet(id='syn25834750')
  if(!is.null(path)){
    setwd(path)
  }
  system(paste('tar -xzvf ',aracne@filePath,sep=''))
    #str1 <- aracne@filePath;
  #str2 <- strsplit(str1,'ARACNE.src.tar.gz')[[1]]
  str3 <- 'ARACNE/'
  cat(str3,'\n')  
  setwd(str3)
  #system(paste(str2,'/ARACNE.src/ARACNE/',sep=''))
  system('make')
  #setwd('../')
  setwd(cwd)
  #return(paste(getwd(),'/ARACNE/',sep=''))
}
