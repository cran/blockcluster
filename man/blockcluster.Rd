\name{blockcluster}
\alias{blockcluster}
\title{Co-Clustering Package}
\description{
  This package performs co-clustering of Binary,
  contingency and categorical datasets.
}
\details{
  This package performs co-clustering of Binary,
  contingency and categorical datasets with utility
  functions to visualize the co-clustered data. The package
  contains a function \code{\link{cocluster}} which perform
  coclustering on various kinds of datasets and returns
  object of appropriate class (refer to documentation of
  \code{\link{cocluster}}). The package also contains
  function \code{\link{cocluststrategy}} (see documentation
  of function to know various slots) which returns an
  object of class \code{\linkS4class{strategy}}. This
  object can be given as input to \code{\link{cocluster}}
  to control various co-clustering paramteres.

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

# A little advanced example with simulated gaussian data
data(gaussiandata)
#set strategy , see documentation for 'cocluststrategy' function for more details.
newstrategy<-cocluststrategy(nbxem=5,nbtry=2,algo="XEMStrategy")
# calling cocluster function with newstrategy
out<-cocluster(gaussiandata,datatype="continuous",nbcocluster=c(2,3),strategy=newstrategy)

# A simple example with simulated contingency data
data(contingencydataunknown)
out<-cocluster(contingencydataunknown,datatype="contingency",nbcocluster=c(2,3))
summary(out)
plot(out)

#Default strategy values
strategy<-cocluststrategy()
summary(strategy)
}

