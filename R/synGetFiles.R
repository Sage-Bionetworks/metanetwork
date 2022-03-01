#' This function pulls files from a synapse project
#' 
#' This function downloads all of the network files in a synapse parent folder using
#' the synID of the folder containing folders of specific network types. It uses
#' the name of the folder plus a pattern id eg. 'Network' to find network files. 
#' for example the `ridgeAIC/` folder in `project_id` will need to contain the
#' `ridgeAICNetwork` file if pattern_id = 'Network' and `project_id` is set to 
#' the parent folder containing `ridgeAIC/`. The files will be downloaded to the
#' path specified in `downloadLocation`.
#' 
#' @param project_id Required. A character vector of a synapse ID of a synapse project.
#' @param pattern_id Required. A character vector too match in the file names of 
#' the children entities in `project_id`.
#' @param downloadLocation Optional. Local directory to download the files to. 
#' Default = `getwd()`
#' 
#' @export 
#' @return A character vector of file paths
#' 
synGetFiles <- function(project_id, pattern_id, downloadLocation = getwd()){
  
  child_obj <- synapser::synGetChildren(project_id, includeTypes = list('folder'))              
  child_list <- as.list(child_obj)
  child_names = c()
  for(k in 1:length(child_list)){
    temp_name = child_list[k]
    temp_name = unlist(temp_name)
    nn = temp_name['name']
    nn = as.character(nn)
    child_names = c(child_names,nn)
  }
  #child_names <- lapply(child_list, `[[`, 1)
  #child_names <- unlist(child_names)
  out_list <- list()
  for (ent in 1:length(child_names)){
    temp_l = synapser::synGetChildren(child_list[[ent]]$id)
    temp_list = temp_l$asList()
    temp_names = c()
    temp_ids = c()
    for(k in 1:length(temp_list)){
      temp_name = temp_list[k]
      temp_name = unlist(temp_name)
      nn = temp_name['name']
      nn = as.character(nn)
      temp_names = c(temp_names,nn)
      nn = as.character(temp_name['id'])
      temp_ids = c(temp_ids,nn)
    }
    # temp_names <- lapply(temp_list, `[[`, 1)
    # temp_names <- unlist(temp_names)
    temp_name_search = paste0(child_names[[ent]],pattern_id)
    temp_names_ids <- grep(temp_name_search, temp_names)

    cat(paste0(temp_names[temp_names_ids],' \n'))
    id_t = temp_ids[temp_names_ids]
    if(length(id_t) == 0){
      cat(paste0('No match for ',temp_name_search,' \n'))
    }else{
      #cat(paste0(id_t,' \n'))
      temp <- synGet(id_t, downloadLocation = downloadLocation)
      out_list <- append(out_list, temp)
    }
  }
  
  if(length(out_list== 0)){
    print("Check your project ID and pattern for input")
  } else{
    print("Downloaded all required network files")
  }
  return(out_list)
}
