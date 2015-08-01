rownameToFirstColumn <- function(DF,colname){
  DF <- as.data.frame(DF)
  DF[,colname] <- row.names(DF)
  DF <- DF[,c(dim(DF)[2],1:(dim(DF)[2]-1))]
  return(DF)
}
