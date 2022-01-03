synGetFiles <- function(project_id, pattern_id, downloadLocation = getwd()){
  
  child_obj <- synapser::synGetChildren(project_id)              
  child_list <- child_obj$asList()
  child_names <- lapply(child_list, `[[`, 1)
  child_names <- unlist(child_names)
  child_names_ids <- grep(pattern_id, child_names)
  out_list <- c()
  # [1] "Filtered Spearman Correlation Table"
  for( ent in child_names_ids){
    message(child_list[[ent]]$name)
    temp <- synGet(child_list[[ent]]$id, downloadLocation =downloadLocation)
    out_list <- append(out_list, temp)
  }
  if(is.na(out_list)){
    print("Check your project ID and pattern for input")
    } else{
  print("Downloaded all required network files")
    }
  return(out_list)
}
