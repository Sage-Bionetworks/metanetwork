require(synapseClient)
synapseLogin()
source('http://depot.sagebase.org/CRAN.R')
pkgInstall("synapseClient")

SVAGeneExpression <- synGet('syn2757147')
cmcSVAGeneExpression <- read.delim(SVAGeneExpression@filePath)

SVATFExpression <- synGet('syn2757149')
cmcSVATFExpression <- read.delim(SVATFExpression@filePath)

data <- rbind(cmcSVAGeneExpression,cmcSVATFExpression)
data <- t(data)
data <- scale(data)

i <- 1000
