#!/usr/bin/env Rscript

# Function to calculate node and network properties (from synapse as RData file)
# Get arguments from comman line
args = commandArgs(TRUE)

# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/home/ec2-user/Work/Github/metanetwork/R')
############################################################################################################

############################################################################################################
#### Libraries ####
library(synapseClient)
library(dplyr)
library(WGCNA)
library(tools)
library(stringr)
library(igraph)
library(data.table)
library(biomaRt)

# Needs the dev branch
library(rGithubClient)
setGithubToken('445838fcb74d2ff824be0cee7c52567000b5f39a')

# login to synapse
synapseLogin()
############################################################################################################

############################################################################################################
#### Github commit ####
# Get github links for provenance
thisFileName = 'nodeNetworkProperties.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/metanetwork", 
                    ref="branch", 
                    refName='nodeProps')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('R/', thisFileName))

# Synapse specific parameters
activityName = 'Node and network properties'
activityDescription = 'Calculating node and network properties of the network'
############################################################################################################

############################################################################################################
#### Get network ####
# Get network from synapse (rda format)
NET_OBJ = synapseClient::synGet(args[1])
FNAME = tools::file_path_sans_ext(NET_OBJ$properties$name)
parentId = NET_OBJ$properties$parentId

# Load sparse network
load(NET_OBJ@filePath)

# Convert lsparseNetwork to igraph graph object
g = igraph::graph.adjacency(sparseNetwork, mode = 'undirected', weighted = NULL, diag = F)
############################################################################################################

############################################################################################################
#### Calculate node and network properties ####
Node.properties <- data.frame(row.names=V(g)$name)

## Basic properties
try({
  # In degree
  writeLines('Calculating in degree...')
  Node.properties$inDegree <- degree(g,
                                     v=V(g), 
                                     mode = 'in', 
                                     loop = T,
                                     normalized = T)
  
  # Out degree
  writeLines('Calculating out degree...')
  Node.properties$outDegree <- degree(g,
                                      v=V(g), 
                                      mode = 'out', 
                                      loop = T,
                                      normalized = T)
  
  # Total degree: An important node is involved in a large number of interactions
  writeLines('Calculating total degree...')
  Node.properties$totalDegree <- degree(g,
                                        v=V(g), 
                                        mode = 'all', 
                                        loop = T,
                                        normalized = T)
  
  # Density: Shows how sparse or dense a network is, for full network the value is 1 and for empty network the value is NaN(0/0)
  writeLines('Calculating density...')
  density <- graph.density(g)
  
  # Diameter: Maximum distance between all the vertex in network
  writeLines('Calculating diameter')
  diameter <- diameter(g)
  
  # Average node degree: Average degree of all nodes in network
  writeLines('Calculating average node degree...')
  avgNodeDegree <- mean(Node.properties$totalDegree, na.rm=T)
  
  # Page Rank:  An important node is likely to receive more connections from other nodes.
  writeLines('Calculating page rank...')
  Node.properties$pageRank <- page.rank(g,
                                        vids = V(g))$vector
  
  # Centralisation: Measures how well network has a star-like topology (closer to 1 more likely the topology is star-like)
  writeLines('Calculating centralisation...')
  centralisation <- vcount(g)/(vcount(g) - 2) * (max(degree(g))/(vcount(g) - 1) - graph.density(g))
  }, 
  silent = T)

# Path level properties 
try({
  # Betwenness centrality: An important node will lie on a high proportion of paths between other nodes in the network.
  writeLines('Calculating betweenness...')
  Node.properties$betweenness <- betweenness(g,
                                             v= V(g),
                                             normalized = T)
  
  # Average path length: Average distance between any two vertex in network
  avgPathLength <- average.path.length(g)
  }, 
  silent = T)
  
# Cluter level properties
try({
  # Clustering coefficeint
  writeLines('Calculating clustering coefficient...')
  Node.properties$clustCoefficient <- transitivity(g,
                                                   v=V(g),
                                                   type = 'barrat')
  
  # Average clustering coefficeint: Gives first hand information on network clusters
  writeLines('Calculating average clustering coefficient...')
  avgClusteringCoefficient <- mean(Node.properties$clustCoefficient, na.rm=T)
  
  # Nearest Neighbor degree
  writeLines('Calculating k nearest neighborhood degree...')
  Node.properties$nearNeighbor <- graph.knn(g,
                                            vids=V(g))$knn
  
  # Closeness centrality: An important node is typically \u201cclose\u201d to, and can communicate quickly with, the other nodes in the network
  writeLines('Calculating closeness centrality ...')
  Node.properties$closeness <- closeness(g,
                                         v= V(g),
                                         mode = 'all',
                                         normalized = T)
  
  # Eigen centrality: An important node is connected to important neighbours.
  writeLines('Calculating eigen centrality...')
  Node.properties$eigen <- alpha.centrality(g,
                                            nodes=V(g))
  }, 
  silent = T)
############################################################################################################

############################################################################################################
#### Write to synapse ####
# Write results to synapse
write.table(Node.properties, paste(FNAME,'nodeProperties.tsv',sep='.'), sep='\t', row.names=F, quote=F)
NODE_OBJ = File(paste(FNAME,'nodeProperties.tsv',sep='.'), name = paste(FNAME,'Node Properties'), parentId = parentId)
annotations(NODE_OBJ) = annotations(NET_OBJ)
NODE_OBJ@annotations$fileType = 'tsv'

# Network level metrics as annotations
try({
  NODE_OBJ@annotations$avgClusteringCoefficient = avgClusteringCoefficient
  NODE_OBJ@annotations$avgNodeDegree = avgNodeDegree
  NODE_OBJ@annotations$avgPathLength = avgPathLength
  NODE_OBJ@annotations$centralisation = centralisation
  NODE_OBJ@annotations$density = density
  NODE_OBJ@annotations$diameter = diameter
},
silent =T)

NODE_OBJ = synStore(NODE_OBJ, 
                    executed = thisFile,
                    used = NET_OBJ,
                    activityName = activityName,
                    activityDescription = activityDescription)
############################################################################################################

############################################################################################################
# Write completed files to synapse
write.table(NODE_OBJ$properties$id, file = 'CompletedNodeIDs.txt', sep='\n', append=T, quote=F, col.names=F, row.names=F)
writeLines(paste('Completed',FNAME,'and stored in',NODE_OBJ$properties$id))
############################################################################################################
