# Function to perform Fishers enrichment analysis
fisherEnrichment <- function(genes.in.significant.set, # A character vector of differentially expressed or some significant genes to test
                             genes.in.gene.set, # A character vector of genes in gene set like GO annotations, pathways etc...
                             genes.in.background # Background genes 
){
  genes.in.significant.set = intersect(genes.in.significant.set, genes.in.background) # back ground filtering
  genes.in.nonsignificant.set = base::setdiff(genes.in.background, genes.in.significant.set)
  genes.in.gene.set = intersect(genes.in.gene.set, genes.in.background) # back ground filtering
  genes.our.gene.set = base::setdiff(genes.in.background,genes.in.gene.set)
  
  tp = length(intersect(genes.in.gene.set, genes.in.significant.set))             
  fp = length(intersect(genes.in.gene.set, genes.in.nonsignificant.set))
  fn = length(intersect(genes.our.gene.set, genes.in.significant.set))
  tn = length(intersect(genes.our.gene.set, genes.in.nonsignificant.set))
  
  pval = fisher.test(matrix(c(tp, fp, fn, tn),nrow=2, ncol=2), alternative="greater")
  odds = (tp*tn)/(fp*fn)
  
  return(data.frame(pval = pval$p.value,
                    ngenes = length(genes.in.gene.set),
                    noverlap = tp,
                    odds = odds,
                    Genes = paste(intersect(genes.in.gene.set, genes.in.significant.set), collapse = '|')
  )
  )
}