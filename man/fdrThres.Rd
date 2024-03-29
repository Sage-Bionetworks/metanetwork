% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fdrThres.R
\name{fdrThres}
\alias{fdrThres}
\title{FDR Threshold}
\usage{
fdrThres(pval, fdr = 0.05)
}
\arguments{
\item{pval}{Required. A vector of uncorected P-Values.}

\item{fdr}{Optional. desired FDR cutoff. (Default = 0.05)
as y}
}
\value{
Corrected PValues
}
\description{
This function applies a user FDR threshold to input p-values
}
