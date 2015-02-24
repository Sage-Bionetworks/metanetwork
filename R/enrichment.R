#function to run simple enrichments
enrichment <- function(list1,list2,list3){
  n2 <- length(list3);
  n1 <- length(list1);
  k <- length(list2);
  a1 <- intersect(list1,list3);
  m <- length(a1);
  n <- n2-m;
  #cat(n2,n1,a1,m,n,'\n')
  
  #s1 <- unlist(lapply(1:n2,getSum,list2,list1));
  s1<-sum(list2%in%list1);
  
  #cat(n2,n1,a1,m,n,s1)
  enr <- (s1/(k))/(m/n2)
  pval <- phyper(q=s1-1,m=m,n=n,k=k,lower.tail=F);
  return(list(enr=enr,pval=pval)); 
}