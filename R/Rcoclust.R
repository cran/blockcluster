#'
#' Co-Clustering Package	
#' 
#' This package performs Co-clustering of Binary, Contingency and Continuous data-sets. 
#' 
#' This package performs Co-clustering of Binary, Contingency and Continuous data-sets with utility functions to
#' visualize the Co-clustered data. The package contains a function \code{\link{cocluster}} which perform Co-clustering
#' on various kinds of data-sets and returns object of appropriate class (refer to documentation of \code{\link{cocluster}}).
#' The package also contains function \code{\link{cocluststrategy}} (see documentation of function to know various slots)
#' which returns an object of class \code{\linkS4class{strategy}}. This object can be given as input to 
#' \code{\link{cocluster}} to control various Co-clustering parameters. Please refer to testmodels.R file which is included in "test" directory
#' to see examples with various models and simulated data-sets.
#' 
#' The package also provide utility functions like summary() and plot() to summarize results and plot
#' the original and Co-clustered data respectively.
#' 
#' 
#' @examples 
#' 
#' # Simple example with simulated binary data
#' #load data
#' data(binarydata)
#' #usage of cocluster function in its most simplest form
#' out<-cocluster(binarydata,datatype="binary",nbcocluster=c(2,3))
#' #Summarize the output results
#' summary(out)
#' #Plot the original and Co-clustered data 
#' plot(out)
#' 
#' 
#' 
#' @name blockcluster
#' @rdname blockcluster
#' 
NULL

#' 
#' An EM strategy to obtain a good optimum.
#' 
#' In Co-clustering, there could be many local optimal where the algorithm may get struck resulting in sub-optimum
#' results. Hence we applied a strategy called XEM strategy to run the EM algorithm. The various steps are defined
#' as follows:
#' 
#' \describe{
#' \item{Step-1, "xem" step:}{Do several runs of: "initialization followed by short run of algorithm (few iterations/high tolerance)". This 
#' parameter is named as "nbxem" in \code{\link{cocluststrategy}} function. Default value is 5. We call this step as xem step.}
#' \item{Step-2, "XEM" step:}{Select the best result of step 1 and make long run of Algorithm(high iterations/low tolerance).We call this step as XEM step.}
#' \item{Step-3}{Repeat step 1 and 2 several times and select the best result. The number of repetitions can be modified via parameter "nbtry"
#' of \code{\link{cocluststrategy}} function. Default value is 2. } 
#' }
#' 
#' 
#' @name XEMStrategy
#' @rdname XEMStrategy
#' 
NULL