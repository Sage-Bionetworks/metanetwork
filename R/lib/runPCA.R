# Function to run principal component analysis
runPCA <- function(genesBySamples, SCALE_DATA_FOR_PCA = TRUE, MIN_PVE_PCT_PC = 1.0) {
  
  # estimate variance in data by PC:
  pca.res <- prcomp(t(genesBySamples), center=TRUE, scale=SCALE_DATA_FOR_PCA, retx=TRUE)
  
  # examine how much variance is explained by PCs, and consider those with PVE >= (MIN_PVE_PCT_PC %):
  pc.var <- pca.res$sdev^2
  pve <- 100 * (pc.var / sum(pc.var))  
  npca <- max(1,length(which(pve >= MIN_PVE_PCT_PC)))
  
  samplePCvals <- pca.res$x[, 1:npca, drop=FALSE]
  
  list(samplePCvals=samplePCvals, pve=pve[1:npca])
}
