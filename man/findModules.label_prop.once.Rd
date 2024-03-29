% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findModules.label_prop.once.R
\name{findModules.label_prop.once}
\alias{findModules.label_prop.once}
\title{Find Modules with Graph Adjacency}
\usage{
findModules.label_prop.once(adj, min.module.size)
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
Deploys `igraph::graph.adjacency` to find module structure.
}
