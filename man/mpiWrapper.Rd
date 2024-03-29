% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mpiWrapper.R
\name{mpiWrapper}
\alias{mpiWrapper}
\title{Runs Sparrow Regression}
\usage{
mpiWrapper(
  data,
  nodes,
  pathv,
  regressionFunction,
  outputpath,
  eigen = NULL,
  regulatorIndex = NULL,
  hosts = NULL
)
}
\arguments{
\item{data}{Required. Expression matrix to be used for network construction.}

\item{nodes}{Required. The number of nodes to run the process over.}

\item{pathv}{Required.}

\item{regressionFunction}{Required. Character vector specifying the specific 
sparrow regression function to deploy. Options are one of c("sparrowZ", 
"sparrow2Z", "sparrow2ZFDR")}

\item{outputpath}{Required. The output path to save the resulting coexpression
network}

\item{eigen}{Optional. (Default = NULL)}

\item{regulatorIndex}{Optional. (Default = NULL)}

\item{hosts}{Optional. (Default = NULL)}
}
\value{
NULL. Saves a sparrow network object to paste0(`outputpath`,
`regressionFunction`,'Network.csv')
}
\description{
This function Sparrow Regression wrapped with the Rmpi::mpi parralel implementation.
User specifies between sparrowZ, sparrow2Z, or sparrow2ZFDR implementation.
}
