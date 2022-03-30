library(dplyr, quietly = TRUE)
library(glmnet, quietly = TRUE)
library(randomForest, quietly = TRUE)
library(Hmisc, quietly = TRUE)
library(lars, quietly = TRUE)
library(WGCNA, quietly = TRUE)
library(synapser, quietly = TRUE)
library(metanetwork, quietly = TRUE)
library(githubr, quietly = TRUE)
library(c3net, quietly = TRUE)
library(config, quietly = TRUE)
library(optparse, quietly = TRUE)
library(data.table, quietly = TRUE)
library(parmigene, quietly = TRUE)
library(reader, quietly = TRUE)
library(Rmpi)
library(parallel)
library(doParallel)
#library(utilityFunctions) -->installation error

# Obtaining the data - From User --------------------------------------------

option_list <- list(make_option(c("-u","--synapse_user"), type="character", action = "store",
                                help = "Synapse User name"),
                    make_option(c("-p","--synapse_pass"), type="character", action = "store",
                                help = "Synapse User Password"),
                    make_option(c("-c","--config_file"), type="character", action = "store",
                                help = "Path to the complete config file"))          
req_args <- parse_args(OptionParser(option_list=option_list))

# Obtaining the data - From Synapse --------------------------------------------

#Setting up the cofig file 
Sys.setenv(R_CONFIG_ACTIVE = "default")
config <- config::get(file = req_args$config_file)
setwd(config$input_profile$temp_storage_loc)

#Linking with Project
synLogin(email = req_args$synapse_user, password = req_args$synapse_pass)
project = Project(config$input_profile$project_id)
project <- synStore(project)

# Data
synID_input = config$input_profile$input_synid
data = synGet(synID_input, downloadLocation = config$input_profile$temp_storage_loc)

# Registering the parallel clusters
if(config$computing_specs$medium_ncores>0){
  nslaves = config$computing_specs$medium_ncores
  mpi.spawn.Rslaves(nslaves=nslaves,hosts=NULL);
}
if(config$computing_specs$heavy_ncores>0){
  nslaves = config$computing_specs$heavy_ncores
  mpi.spawn.Rslaves(nslaves=nslaves,hosts=NULL);
}

# Performing the analysis -------------------------------------------------

net_methods = config$input_profile$network_method
data = reader::reader(data$path)

if( is.null(config$input_profile$na_fill)){
  print('Data not normalized for missing values. Ignore if using mrnet method.')
}else{
   if (config$input_profile$na_fill == 'Winsorize'){
    data <- metanetwork::winsorizeData(data)
  }
}

for (method in net_methods){# Assuming we have more methods - not developing for now

  switch(method,
         "c3net" = c3netWrapper(data, pval = config$input_profile$p_val_c3net, config$output_profile$output_path),# What does this path define in main function?
         "mrnet" = mrnetWrapper(data=data, path = NULL, pval = config$input_profile$p_val_mrnet, outputpath=config$output_profile$output_path,  tool_storage_loc = config$input_profile$temp_storage_loc),
         "wgcna" = wgcnaTOM(data=data, path = NULL, pval = config$input_profile$p_val_wgcna, outputpath=config$output_profile$output_path, 
                            config$input_profile$rsquaredcut, config$input_profile$defaultnaPower),
         
         "lassoAIC" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                 outputpath = config$output_profile$output_path),
         "lassoBIC" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                 outputpath = config$output_profile$output_path),
         "lassoCV1se" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                   outputpath = config$output_profile$output_path),
         "lassoCVmin" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                   outputpath = config$output_profile$output_path),
         "ridgeAIC" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                 outputpath = config$output_profile$output_path),
         "ridgeBIC" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                 outputpath = config$output_profile$output_path),
         "ridgeCV1se" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                   outputpath = config$output_profile$output_path),
         "ridgeCVmin" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                   outputpath = config$output_profile$output_path),
         "sparrowZ" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                 outputpath = config$output_profile$output_path),
         "sparrow2Z" = mpiWrapper(data, nodes = config$computing_specs$medium_ncores, pathv = NULL, regressionFunction = method,
                                  outputpath = config$output_profile$output_path),
         "genie3" = mpiWrapper(data, nodes = config$computing_specs$heavy_ncores, pathv = NULL, regressionFunction = method,
                               outputpath = config$output_profile$output_path),
         "tigress" = mpiWrapper(data, nodes = config$computing_specs$heavy_ncores, pathv = NULL, regressionFunction = method,
                                outputpath = config$output_profile$output_path))
  
  output_filename <- list.files(pattern = method )
}
if(config$computing_specs$heavy_ncores>0){
  mpi.close.Rslaves()
}

if(config$computing_specs$medium_ncores>0){
  mpi.close.Rslaves()
}
  

# Obtaining the data - For provenance --------------------------------------------

all.annotations <- synGetAnnotations(config$input_profile$input_synid)

checkAnnotations <- function(annotations, config){
  annot_default <- list(
    dataType = NULL,
    resourceType = NULL,
    metadataType = NULL,
    isModelSystem = NULL,
    isMultiSpecimen = NULL,
    fileFormat = NULL,
    grant = NULL,
    species = NULL,
    organ = NULL,
    tissue = NULL,
    study = NULL, 
    consortium = NULL,
    assay = NULL
  )
  for (item in names(annot_default)){
    if (!is.null(config$provenance$annotations$item)){
      annot_default$item = config$provenance$annotations$item[[1]]
    }
    else if (!is.null(annotations$item)){
      annot_default$item = annotations$item[[1]]
    }
  }
}

all.annotations <- checkAnnotations(all.annotations,config)

# Pull Git Provenance
thisRepo = NULL
thisFile = NULL

try(
  thisRepo <- githubr::getRepo(
    repository = config$provenance$code_annotations$repository,
    ref = config$provenance$code_annotations$ref,
    refName = config$provenance$code_annotations$ref_name
  ), silent = TRUE
)
try(
  thisFile <- githubr::getPermlink(
    repository = thisRepo,
    repositoryPath = config$provenance$code_annotations$repository_path
  ), silent = TRUE
)

## Store the file
activity <- synapser::synGet(config$input_profile$project_id)

dataFolder <- Folder(method,parent = config$input_profile$project_id)
dataFolder <- synStore(dataFolder)

#### - push the config file req_args$config_file
#split path to find file name

Config_OBJ <- synapser::synStore( synapser::File( 
    path = req_args$config_file,
    name = tail(strsplit(req_args$config_file,'/')[[1]], n =1),
    parentId = dataFolder$properties$id),
    used = config$input_profile$input_synid,
    activityName = config$provenance$activity_name,
    executed = thisFile,
    activityDescription = config$provenance$activity_description
  )

syn_config <- Config_OBJ$properties$id

####

if( method == 'wgcna') {
  filePath <- gsub(
    '//',
    '/', 
    paste0( 
      config$output_profile$output_path,
      '/',
      list.files(config$output_profile$output_path, pattern = method ))
    )
  mdpath <- gsub('.txt', 'md5.out', filePath)
  mdpath <- gsub('.csv', 'md5.out', mdpath)
  syn_name <- gsub(
    '\\.txt',
    '', 
    gsub(
      '\\.csv',
      '',
      list.files(config$output_profile$output_path, pattern = method ))
    )
  syn_name <- gsub('wgcna', 'wgcna ', gsub( 'Network',' Network',syn_name))
  syn_name <- gsub('Power', 'Power ', gsub( 'Soft','Soft ',syn_name))
  syn_name <- gsub('Overlap', ' Overlap ', syn_name)
}else{
  filePath <- gsub(
    '//',
    '/', 
    paste0(
      config$output_profile$output_path, 
      '/', 
      config$input_profile$network_method,'Network.csv')
  )
  mdpath <- config$output_profile$md5_output_path
  syn_name <- config$output_profile$output_name
}

for(file in filePath) {
  ENRICH_OBJ <- synapser::synStore( synapser::File( 
    path = file,
    name = syn_name[which(filePath == file)],
    parentId = dataFolder$properties$id),
    used = c( config$input_profile$input_synid, syn_config),
    activityName = config$provenance$activity_name,
    executed = thisFile,
    activityDescription = config$provenance$activity_description
  )

  # Formatting the network md5 ------------------------------------------------------

  md5Command <- paste0('md5sum ', file)
  md5 <- strsplit(system(md5Command, intern = TRUE), '  ')[[1]][1]
  cat(md5, '\n', file = mdpath[which(filePath == file)], sep = '')

  ENRICH_OBJ <- synapser::synStore( synapser::File( 
    path = mdpath[which(filePath == file)],
    name = gsub(
      config$output_profile$output_path,
      '',
      gsub('\\.out','',mdpath[which(filePath == file)])),
    parentId = dataFolder$properties$id),
    used = c( config$input_profile$input_synid, syn_config),
    activityName = config$provenance$activity_name,
    executed = thisFile,
    activityDescription = config$provenance$activity_description
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
}

if((!is.na(config$computing_specs$heavy_ncores)) || (!is.na(config$computing_specs$medium_ncores))){
  Rmpi::mpi.quit(save = "no")
}

