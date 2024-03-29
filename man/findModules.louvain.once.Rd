% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findModules.louvain.once.R
\name{findModules.louvain.once}
\alias{findModules.louvain.once}
\title{Find Modules with Network Adjacency Matrix Using Louvain Clustering}
\usage{
findModules.louvain.once(adj, min.module.size)
}
\arguments{
\item{adj}{A n x n upper triangular adjacency in the matrix class format.}

\item{min.module.size}{Optional. Integer between 1 and n genes. (Default = 30)}
}
\value{
GeneModules = n x 3 dimensional data frame with column names as Gene.ID,
moduleNumber, and moduleLabel.
}
\description{
This function tries to get modules from network adjacency matrix using igraph's 
louvain clusting function.
}
