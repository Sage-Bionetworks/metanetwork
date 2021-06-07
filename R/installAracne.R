installAracne <- function(path=NULL){
  arc_path = system.file("tools", "ARACNE.src.tar.gz", package = "metanetwork")
  system(paste('tar -xzvf ',arc_path,sep=''))
  #str1 <- aracne@filePath;
  #str2 <- strsplit(str1,'ARACNE.src.tar.gz')[[1]]
  temp_path = paste0(config$input_profile$temp_storage_loc,"/ARACNE")
  setwd(temp_path)
  cat(temp_path,'\n')  
  setwd(temp_path)
  #system(paste(str2,'/ARACNE.src/ARACNE/',sep=''))
  system('make')
  #setwd('../')
  #return(paste(getwd(),'/ARACNE/',sep=''))
}
