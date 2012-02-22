\docType{methods}
\name{plot}
\alias{plot}
\alias{plot,BinaryOptions-method}
\alias{plot,ContingencyOptions-method}
\alias{plot,ContinuousOptions-method}
\title{Plot function.}
\arguments{
  \item{x}{output object from \code{\link{cocluster}}.}

  \item{y}{Ignored}

  \item{...}{Additional argument(s) . Currently we support
  only one additional argument which is "asp". If this is
  set to FALSE the original aspect ratio is not conserved.
  In this case, R graphics will display the data in a way
  that fits the screen and is suitable for viewing. By
  default "asp" is TRUE.}
}
\description{
  This funtion plot the original and co-clustered datasets.
}

