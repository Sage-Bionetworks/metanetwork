#install aracne function
installAracne <- function(path=NULL){

  cwd <- getwd()
    #str1 <- aracne@filePath;
  #str2 <- strsplit(str1,'ARACNE.src.tar.gz')[[1]]
    str3 <- '~/ARACNE/'
    cat(str3,'\n')  
    setwd(str3)
    #system(paste(str2,'/ARACNE.src/ARACNE/',sep=''))
    system('make')
    #setwd('../')
    setwd(cwd)
    #return(paste(getwd(),'/ARACNE/',sep=''))
   
}
