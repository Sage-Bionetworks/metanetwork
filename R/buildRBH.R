library(dplyr)

if(!("fastICA" %in% rownames(installed.packages()))){
  install.packages("fastICA")
}
library(fastICA)

synapser::synLogin()
whole <- read.csv(synapser::synGet('syn21266454')$path, row.names = 1)
whole <- read.csv('~/Documents/Test_Folder/C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv')

#Randomize the columnsa
set.seed(42)
whole_rand <- whole[ ,
                     sample(colnames(whole),
                            size = dim(whole)[2],
                            replace = FALSE
                           )
]

T_1 <- whole_rand[,
                  1:floor(ncol(whole_rand)/2)
                 ]
T_2 <- whole_rand[,
                  (floor(ncol(whole_rand)/2)+1):ncol(whole_rand)
                  ]

T_1[is.na(T_1)] = 0
T_2[is.na(T_2)] = 0

#MF
mf_1 <- fastICA::fastICA(T_1, 100)
mf_2 <- fastICA::fastICA(T_2,100)

Metagenes_1 <- mf_1$S
Metagenes_2 <- mf_2$S
MetaSamples_1 <- mf_1$A
MetaSamples_2 <- mf_2$A

# These inner Loop could a stand alone function implemented with parApply
# Returns a list of each run, and convert to df with do.call(rbind,<list>)
#' @param i A Vector to compare to each vector in j
#' @param j A data frame of vectors to compare to i
#' @param method  a correlation method ie "pearson", "kendall", or "spearman"
#' @param name_ind a numerical index id for column names ie. 1 or 2. Default = 1
cor_func <- function( i, j, method, name_ind = 1) {
  old_cor = 0
  idx = 0
  p_vals <- NULL
  cors <- NULL
  res <- c(NA, NA, NA ,NA)
  
  for (comp_row in 1:nrow(j)){
    temp_cor = cor(i, j[comp_row,], method = method)
    cors <- c(cors,temp_cor)
    
    temp_cortest = cor.test(i, j[comp_row,], method = method)
    p_vals <- c(p_vals, temp_cortest$p.value)
    #_# Need to take Mod of both vals to ensure absolute value is taken 
    if (Mod(temp_cor) > Mod(old_cor) ){
      old_cor = temp_cor
      idx = comp_row
    }
  }
  
  adj_p <- p.adjust(p_vals, method = 'fdr' )
  
  if( length(cors[ which(adj_p < 0.05) ]) > 0 ){
    keep_ind <- which(abs(cors) == max(abs(cors[ which(adj_p < 0.05) ])))
    keep <- cors[keep_ind]
  }else{
    keep <- NA
    keep_ind <- NA
  }
  
  res[1] = idx
  res[2] = old_cor
  res[3] = keep_ind
  res[4] = keep
  names(res) <- c(paste0('Index_',name_ind),
                  paste0('Correlation_',name_ind),
                  paste0('Index_adj_',name_ind),
                  paste0('Correlation_adj_',name_ind)
  )
  return(res)
}

cl <- parallel::makeCluster(parallel::detectCores()-1)

start <- Sys.time()
res <-as.data.frame(t(as.matrix(parallel::parApply(
  cl = cl, Metagenes_1, 1, cor_func,
  j = Metagenes_2, method = "pearson", name_ind = 1
))))
Sys.time()-start

start2 <- Sys.time()
# Run inverse:
res_2 <- as.data.frame(t(as.matrix( parallel::parApply(
  cl = cl, Metagenes_2, 1, cor_func,
  j = Metagenes_1, method = "pearson", name_ind = 2
))))
Sys.time()-start2


## I Get lost here: 
#res_2_mod = res_2[res$Index_1,]
res_final <- cbind(res,res_2)
sink_final <- res_final

#idx_rbh = which(res_final$Genes == res_final$Index_2)
idx_rbh = which(res_final$Index_1 == res_final$Index_2)

res_rbh = res_final[idx_rbh,]
#                 Index_1 Correlation_1 Index_2 Correlation_2
# PDE1B|Q01064        363    0.42195208     363    0.45691483
# HEXA|P06865        8817    0.14453328    8817   -0.15883453
# TMEM63C|Q9P1W3      363    0.52228077     363    0.51437329
# ABHD8|Q96I13       8817    0.05199916    8817    0.04799616
# EPM2AIP1|Q7L775    4035    0.43378938    4035    0.46993216
# L1CAM|P32004        363    0.43625009     363    0.39460592
# PPP1R12B|O60237    8817    0.06863327    8817    0.01856145
# AMPH|P49418         363    0.43281238     363    0.40487099
# PHAX|Q9H814        6742    0.42319546    6742    0.37327267
# LASP1|Q14847       8817   -0.15735168    8817    0.16473248
# FADS2|O95864       8583    0.28780491    8583    0.38937385
# TP53I3|Q53FA7      8817   -0.11823845    8817    0.13080624
# CAMK2N1|Q7Z7J9     8817    0.04577299    8817    0.12704475
# PCSK2|P16519        363    0.54782769     363    0.57389052

res_final[ order(abs(res_final$Correlation_1), decreasing = T),][1:20,]
#                    Index_1 Correlation_1 Index_2 Correlation_2
# CAP1|Q01518           2795     0.9378616    7153     0.4333556
# MT1F|P04733           7671     0.8327836    1571     0.5173576
# HBA2|P69905           6256     0.8246218      95     0.6082767
# DHX57|Q6P158           674     0.8228866    6099     0.3976221
# HBB|P68871            6256     0.8116894      95     0.6491741
# DDX41|J3KNN5          2614     0.8097495    5297     0.6701653
# MT1E|P04732           7671     0.8096307    1571     0.5131141
# MT2A|P02795           7671     0.7942879    6770     0.5362846
# HBD|P02042            6256     0.7799367      95     0.6169417
# HBG2|P69891           7569     0.7754258     748     0.8097495
# GPRC5B|Q9NZH0         2506     0.7738509    6781     0.3457700
# FAM219A|A0A0A0MRW3    2614     0.7706290    5297     0.5786435
# SLC4A1|P02730         6256     0.7693038      95     0.6155648
# SMC3|Q9UQE7           2614     0.7676008    5296     0.4520630
# HBG2|P69892           8521     0.7664621    5230     0.7245093
# KRT5|P13647           3187     0.7604107    5641     0.8229051
# GOLGB1|Q14789-2        674     0.7517579    5297     0.4414087
# SMYD2|Q9NRG4          2614     0.7509821    8798     0.3994942
# HLA-A|P30450          1039     0.7488569    8063     0.3013216
# PAIP2|Q9BPZ3          2614     0.7474692    8783     0.3003472

res_final[ order(abs(res_final$Correlation_2), decreasing = T),][1:20,]
#               Index_1 Correlation_1 Index_2 Correlation_2
# KRT1|P04264       1529     0.7202377    5641     0.9378616
# KRT10|P13645      3632     0.6655440    5641     0.8635552
# ASAP1|H0YBM4      2090     0.5807814    3868     0.8327836
# OPTN|Q96CV9-2     3985     0.3945329     204     0.8246218
# KRT5|P13647       3187     0.7604107    5641     0.8229051
# PPP3R1|D3YTA9     6924     0.5255346    1964     0.8228866
# HBG2|P69891       7569     0.7754258     748     0.8097495
# KRT9|P35527       1689     0.6605322    5641     0.7780387
# COL14A1|Q05707    8721     0.2823754    2614     0.7754258
# MFAP2|P55001      8735     0.2801762    2614     0.7753974
# FLNA|P21333       3434     0.6528353    5215     0.7738509
# PKP2|Q99959-2     3313     0.3819015    4863     0.7664621
# KRT14|P02533      1689     0.6965367    5641     0.7608946
# SYT2|Q8N9I0       7052     0.5784543    8620     0.7604107
# TPM2|P07951       4691     0.6683940    5215     0.7603259
# FHOD3|Q2V2M9-4    7394     0.4047139    5297     0.7488569
# SCN4B|Q8IWT1      7052     0.5585262    8620     0.7358570
# TPM2|P07951-3     4691     0.6493001    5215     0.7353042
# EMILIN1|Q9Y6C2    4933     0.3982774    2614     0.7352527
# KRT2|P35908       3632     0.6226017    5641     0.7346409

max_rdx = which(res_rbh$Correlation_1 == max(res_rbh$Correlation_2))
res_rbh[max_rdx,]

gene_names = whole$X
head(res_rbh)
gene_names[5]

result_rbh = res_rbh[,1:3]
colnames(result_rbh) = c("Genes_from_Set1", "Genes_from_Set2","Correlation")
for ( k in 1:nrow(result_rbh)){
  temp = as.integer(result_rbh$Genes_from_Set1[k])
  result_rbh$Genes_from_Set1[k] = gene_names[temp]
  temp = as.integer(result_rbh$Genes_from_Set2[k])
  result_rbh$Genes_from_Set2[k] = gene_names[temp]
}

result_rbh = result_rbh %>% arrange(desc(Correlation))
head(result_rbh)


install.packages("MCL")
library(MCL)

length(unique(result_rbh$Genes_from_Set1))
length(unique(result_rbh$Genes_from_Set2))
dim(result_rbh)

mcl(result_rbh)










