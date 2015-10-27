###Function to run genie3 on the data
genie3 <- function(y,x){
  #function to run genie3 on a single gene
  #IncNodePurity
  require(randomForest)
  rf <- randomForest(x,y,mtry=round(sqrt(ncol(x))),ntree=1000,importance=TRUE)
  im <- importance(rf)[,'IncNodePurity'];
  return(im)
}
