#' This function applies ARACNE on the data
#' 
#' This function installs the ARACNE system Function
#' 
#' @param path Optional. Path is currently not a functional argument. This function
#' automatically installs ARACNe from `inst/tools`. Default = NULL.
#' tool_storage_loc Required. Provides the directory inside docker to 
#' temporarily store the ARACNE files and package 
#'
#' @export
#' @return NULL
#' 
installAracne <- function(path=NULL,tool_storage_loc = config$input_profile$temp_storage_loc){
  arc_path = system.file("inst/tools", "ARACNE.src.tar.gz", package = "metanetwork")
  system(paste('tar -xzvf ',arc_path,sep=''))
  #str1 <- aracne@filePath;
  #str2 <- strsplit(str1,'ARACNE.src.tar.gz')[[1]]
  temp_path = paste0(tool_storage_loc,"/ARACNE")
  setwd(temp_path)
  cat(temp_path,'\n')  
  setwd(temp_path)
  #system(paste(str2,'/ARACNE.src/ARACNE/',sep=''))
  system('make')
  #setwd('../')
  #return(paste(getwd(),'/ARACNE/',sep=''))
}
