% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeDriverDistance.R
\name{computeDriverDistance}
\alias{computeDriverDistance}
\title{Compute Graph Distance}
\usage{
computeDriverDistance(geneSet, graph)
}
\arguments{
\item{geneSet}{Required. A user specified gene set corresponding to vertex names
in `graph`. All the names are not required to be entirely represented in the 
graph object.}

\item{graph}{Required. An igraph graph object consisting of vertices of genes
and edges representing co-expression between genes.}
}
\value{
a list consisting of a mean distance of genes within the largest component 
to each other and the paiwise list of distances.
}
\description{
Finds the finds the maximal (weakly or strongly) connected component of `graph` 
then returns the distances between members in `geneSet` which are found in the 
largest component.
}
