\name{blockcluster}
\alias{blockcluster}
\title{Co-Clustering Package}
\description{
  This package performs co-clustering of Binary,
  Contingency and Continuous datasets.
}
\details{
  This package performs co-clustering of Binary,
  Contingency and Continuous datasets with utility
  functions to visualize the co-clustered data. The package
  contains a function \code{\link{cocluster}} which perform
  coclustering on various kinds of datasets and returns
  object of appropriate class (refer to documentation of
  \code{\link{cocluster}}). The package also contains
  function \code{\link{cocluststrategy}} (see documentation
  of function to know various slots) which returns an
  object of class \code{\linkS4class{strategy}}. This
  object can be given as input to \code{\link{cocluster}}
  to control various co-clustering paramteres. Please refer
  to testmodels.R file which is included in "test"
  directory to see examples with various models and
  simulated data-sets.

  The package also provide utility functions like summary()
  and plot() to summarize results and plot the original and
  co-clustered data respectively.
}
\examples{
# Simple example with simulated binary data
#load data
data(binarydata)
#usage of cocluster function in its most simplest form
out<-cocluster(binarydata,datatype="binary",nbcocluster=c(2,3))
#Summarize the output results
summary(out)
#Plot the original and co-clustered data
plot(out)
}

