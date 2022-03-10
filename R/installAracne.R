#' This function applies ARACNE on the data
#' 
#' This function installs the ARACNE system Function
#' 
#' @param path Optional. Path is currently not a functional argument. This function
#' automatically installs ARACNe from `inst/aracne`. Default = NULL.
#' @param tool_storage_loc Required. Provides the directory inside docker to 
#' temporarily store the ARACNE files and package. (Default = config$input_profile$temp_storage_loc)
#' 
#' @export
#' @return NULL
#' 
installAracne <- function(path=NULL, tool_storage_loc){
  if(is.null(path)) {
    arc_path = system.file("inst/aracne", "ARACNE.src.tar.gz", package = "metanetwork")
  }else {
    arc_path = path
  }
  
  #Unzip Aracne into the temp files directory 
  system(paste('tar -xzvf ',arc_path, ' -C ', tool_storage_loc, sep=''))
  temp_path = gsub( '//', '/', paste0(tool_storage_loc,"/ARACNE"))
  #Change to Aracne directory and install
  setwd(temp_path)
  cat(temp_path,'\n')  
  setwd(temp_path)
  
  system('make')
  
}
