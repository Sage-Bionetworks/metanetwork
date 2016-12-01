## Module test with 100 gene network
# Function to get bicNetworks from synapse and find modules and push results back to synapse

# Load libraries
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

library(igraph)
library(metanetwork)

library(synapseClient)

synapseLogin()

# Synapse parameters
net.id = 'syn7291744'

# Get network from synapse
net.obj = synapseClient::synGet(net.id)
load(net.obj@filePath)

# Get the adjacency matrix
adj <- bicNetworks$network[1:100,1:100]

# Perform module identification
source('./R/findModules.edge_betweenness.R')
mod1 = findModules.edge_betweenness(adj)

source('./R/findModules.fast_greedy.R')
mod2 = findModules.fast_greedy(adj)

source('./R/findModules.infomap.R')
mod3 = findModules.infomap(adj)

source('./R/findModules.label_prop.R')
mod4 = findModules.label_prop(adj)

source('./R/findModules.leading_eigen.R')
mod5 = findModules.leading_eigen(adj)

source('./R/findModules.louvain.R')
mod6 = findModules.louvain(adj)

source('./R/findModules.spinglass.R')
mod7 = findModules.spinglass(adj)

source('./R/findModules.walktrap.R')
mod8 = findModules.walktrap(adj)

source('./R/findModules.hclust.R')
mod9 = findModules.hclust(adj)

source('./R/findModules.GANXiS.R')
mod10 = findModules.GANXiS(adj, '/mnt/Github/metanetwork/R/GANXiS_v3.0.2/')

source('./R/findModules.CFinder.R')
mod11 = findModules.CFinder(adj, '/mnt/Github/metanetwork/CFinder-2.0.6--1448/')

source('./R/compute.Modularity.R')
Q1 = compute.Modularity(adj, mod2, method = 'Newman1')
Q2 = compute.Modularity(adj, mod2, method = 'Newman2')

source('./R/compute.LocalModularity.R')
NQ = compute.LocalModularity(adj, mod2)

source('./R/compute.ModularityDensity.R')
Qds = compute.ModularityDensity(adj, mod2)

source('./R/compute.ModuleQualityMetric.R')
Mod.metrics = compute.ModuleQualityMetric(adj, mod2)
