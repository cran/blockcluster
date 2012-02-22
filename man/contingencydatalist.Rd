\docType{data}
\name{contingencydatalist}
\alias{contingencydatalist}
\title{Simulated Contingency Data-set}
\format{A data list consisting of following data:

\describe{

  \item{\code{data}}{A data matrix consisting of 1000 rows and 100 columns.}

  \item{\code{roweffects}}{A numeric vector of size 1000. Each value represent row effect of corresponding row.}

  \item{\code{columneffects}}{A numeric vector of size 100. Each value represent column effect of corresponding column.}
}}
\description{
  It is a contingency data-set simulated using Poisson distribution. The row and column effects is known for 
  this data-set. It consist of two clusters in rows and three clusters in columns. 
}
\examples{
data(contingencydatalist)
}
\keyword{datasets}

