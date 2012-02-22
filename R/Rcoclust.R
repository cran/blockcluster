#'
#' Co-Clustering Package	
#' 
#' This package performs co-clustering of Binary, contingency and categorical datasets. 
#' 
#' This package performs co-clustering of Binary, contingency and categorical datasets with utility functions to
#' visualize the co-clustered data. The package contains a function \code{\link{cocluster}} which perform coclustering
#' on various kinds of datasets and returns object of appropriate class (refer to documentation of \code{\link{cocluster}}).
#' The package also contains function \code{\link{cocluststrategy}} (see documentation of function to know various slots)
#' which returns an object of class \code{\linkS4class{strategy}}. This object can be given as input to 
#' \code{\link{cocluster}} to control various co-clustering paramteres. 
#' 
#' The package also provide utility functions like summary() and plot() to summarize results and plot
#' the original and co-clustered data respectively.
#' 
#' 
#' @examples 
#' # Simple example with simulated binary data
#' #load data
#' data(binarydata)
#' #usage of cocluster function in its most simplest form
#' out<-cocluster(binarydata,datatype="binary",nbcocluster=c(2,3))
#' #Summarize the output results
#' summary(out)
#' #Plot the original and co-clustered data 
#' plot(out)
#' 
#' # A little advanced example with simulated gaussian data
#' data(gaussiandata)
#' #set strategy , see documentation for 'cocluststrategy' function for more details.
#' newstrategy<-cocluststrategy(nbxem=5,nbtry=2,algo="XEMStrategy")
#' # calling cocluster function with newstrategy
#' out<-cocluster(gaussiandata,datatype="continuous",nbcocluster=c(2,3),strategy=newstrategy)
#' 
#' # A simple example with simulated contingency data
#' data(contingencydataunknown)
#' out<-cocluster(contingencydataunknown,datatype="contingency",nbcocluster=c(2,3))
#' summary(out)
#' plot(out)
#' 
#' #Default strategy values
#' strategy<-cocluststrategy()
#' summary(strategy)
#' 
#' @name blockcluster
#' @rdname blockcluster
#' 
NULL

#' 
#' An EM strategy to obtain a good optimum.
#' 
#' In co-clustering, there could be many local optima where the algorithm may get struct resulting in sub-optimum
#' results. Hence we applied a strategy called XEM strategy to run the EM algorithm. The various steps are defined
#' as follows:
#' 
#' \describe{
#' \item{Step-1, "xem" step:}{Do several runs of: "initialization followed by short run of algorithm (few iterations/high tolerance)". This 
#' parameter is named as "nbxem" in \code{\link{cocluststrategy}} function. Default value is 5. We call this step as xem step.}
#' \item{Step-2, "XEM" step:}{Select the best result of step 1 and make long run of Algorithm(high iterations/low tolerance).We call this step as XEM step.}
#' \item{Step-3}{Repeat step 1 and 2 several times and select the best result. The number of repeatitions can be modified via parameter "nbtry"
#' of \code{\link{cocluststrategy}} function. Default value is 2. } 
#' }
#' 
#' 
#' @name XEMStrategy
#' @rdname XEMStrategy
#' 
NULL