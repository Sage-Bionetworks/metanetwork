makeModulePreservationSubmissions <- function(ref, test, refNet, testNet, refModLabels, testModLabels, refExp, testExp, n= 100){
  
  folderName = paste(getwd(), paste(ref,'as_ref',test,'as_test',sep='_'), sep='/')
  system(paste('mkdir',folderName))
  
  # Track all subission scripts in one shell script
  fp_all = file(paste0(folderName, '/allSubmissions.sh'),'w+')    
  cat('#!/bin/bash',file=fp_all,sep='\n')
  close(fp_all)
  
  # Package actual data and submit them to sge
  netData = list(refNet = refNet, testNet = testNet, 
                 refModLabels = refModLabels, testModLabels = testModLabels, 
                 refExp = refExp, testExp = testExp)                                    
  save(list = 'netData', file = paste(folderName, 'Main.RData',sep='/'))
  
  # Create main submission script
  fp = file (paste0(folderName,'/Main.sh'), "w+")
  cat('#!/bin/bash',
      'sleep 30',
      paste('Rscript','/home/ec2-user/Work/Github/metanetwork/R/modulePreservationAnalysis.SGE.R','Main.RData',folderName,'Main'),
      file = fp,
      sep = '\n')
  close(fp)
    
  # Add submission script to allSubmission list
  fp_all = file(paste0(folderName, '/allSubmissions.sh'),'a+')
  cat(paste('qsub','-cwd','-V',paste(folderName,'Main.sh',sep='/'),
            '-o',paste(folderName,'Main.o',sep='/'),
            '-e',paste(folderName,'Main.e',sep='/'),
            '-l mem=7GB'),
      file=fp_all,
      sep='\n')
  close(fp_all)
  
  # Create random networks for sge submission
  for (i in 1:n){
    ind = sample(vcount(refNet))
    adjRefNet = as_adj(refNet)
    refModLabels = refModLabels[rownames(adjRefNet)]
    rownames(adjRefNet)[ind] = rownames(adjRefNet)
    colnames(adjRefNet)[ind] = colnames(adjRefNet)
    permRefNet = igraph::graph.adjacency(adjRefNet, mode = 'undirected', weighted = NULL, diag = F)
    
    permModLabels = refModLabels
    names(permModLabels)[ind] = names(refModLabels)
    
    ind = sample(vcount(testNet))
    adjTestNet = as_adj(testNet)
    testModLabels = testModLabels[rownames(adjTestNet)]
    rownames(adjTestNet)[ind] = rownames(adjTestNet)
    colnames(adjTestNet)[ind] = colnames(adjTestNet)
    permTestNet = igraph::graph.adjacency(adjTestNet, mode = 'undirected', weighted = NULL, diag = F)
    
    # Package actual data and submit them to sge
    netData = list(refNet = permRefNet, testNet = permTestNet, 
                   refModLabels = permModLabels, testModLabels = testModLabels, 
                   refExp = refExp, testExp = testExp)                                    
    save(list = 'netData', file = paste(folderName, paste('Rand',i,'RData',sep='.'),sep='/'))
    
    # Create main submission script
    fp = file (paste(folderName, paste('Rand',i,'sh',sep='.'),sep='/'), "w+")
    cat('#!/bin/bash',
        'sleep 30',
        paste('Rscript','/home/ec2-user/Work/Github/metanetwork/R/modulePreservationAnalysis.SGE.R',paste(folderName, paste('Rand',i,'RData',sep='.'),sep='/'),folderName,paste('Rand',i,sep='.')),
        file = fp,
        sep = '\n')
    close(fp)
    
    # Add submission script to allSubmission list
    fp_all = file(paste0(folderName, '/allSubmissions.sh'),'a+')
    cat(paste('qsub','-cwd','-V',paste(folderName, paste('Rand',i,'sh',sep='.'),sep='/'),
              '-o',paste(folderName, paste('Rand',i,'o',sep='.'),sep='/'),
              '-e',paste(folderName, paste('Rand',i,'e',sep='.'),sep='/'),
              '-l mem=7GB'),
        file=fp_all,
        sep='\n')
    close(fp_all)
  }
}
