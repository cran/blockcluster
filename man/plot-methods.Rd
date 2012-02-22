\docType{methods}
\name{plot}
\alias{plot}
\alias{plot,BinaryOptions-method}
\alias{plot,CategoricalOptions-method}
\alias{plot,ContingencyOptions-method}
\alias{plot,ContinuousOptions-method}
\title{Plot function.}
\arguments{
  \item{x}{output object from \code{\link{cocluster}}.}

  \item{y}{Ignored}

  \item{...}{Additional argument(s) . Currently we support
  two additional argument. "asp": If this is set to TRUE
  the original aspect ratio is conserved. By default "asp"
  is FALSE. "type" : This is the type of plot which is
  either "cocluster" or "distribution". The corresponding
  plots are Co-clustered data and distributions and mixture
  densities for Co-clusters respectively. Default is
  "cocluster" plot.}
}
\description{
  This function plot the original and Co-clustered
  data-sets.
}

