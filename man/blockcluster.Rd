\name{blockcluster}
\alias{blockcluster}
\title{Co-Clustering Package}
\description{
  This package performs Co-clustering of Binary,
  Contingency and Continuous data-sets.
}
\details{
  This package performs Co-clustering of Binary,
  Contingency and Continuous data-sets with utility
  functions to visualize the Co-clustered data. The package
  contains a function \code{\link{cocluster}} which perform
  Co-clustering on various kinds of data-sets and returns
  object of appropriate class (refer to documentation of
  \code{\link{cocluster}}). The package also contains
  function \code{\link{cocluststrategy}} (see documentation
  of function to know various slots) which returns an
  object of class \code{\linkS4class{strategy}}. This
  object can be given as input to \code{\link{cocluster}}
  to control various Co-clustering parameters. Please refer
  to testmodels.R file which is included in "test"
  directory to see examples with various models and
  simulated data-sets.

  The package also provide utility functions like summary()
  and plot() to summarize results and plot the original and
  Co-clustered data respectively.
}
\examples{
# Simple example with simulated binary data
#load data
data(binarydata)
#usage of cocluster function in its most simplest form
out<-cocluster(binarydata,datatype="binary",nbcocluster=c(2,3))
#Summarize the output results
summary(out)
#Plot the original and Co-clustered data
plot(out)
}

