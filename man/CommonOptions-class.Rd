% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optionclasses.R
\docType{class}
\name{CommonOptions-class}
\alias{CommonOptions-class}
\title{Common Input/Output options.}
\description{
This class contains all the input options and common output options for all
kinds of data-sets (Binary, Categorical, Contingency and Continuous).
}
\details{
\describe{
\item{The following are the various input options:}{}
\item{data: }{Input data.}
\item{datatype: }{Type of data.}
\item{semisupervised: }{Boolean value specifying if Co-clustering is semi-supervised or not.}
\item{model: }{Model to be run for co-clustering.}
\item{nbcocluster: }{Number of row and column clusters.}
\item{strategy: }{Input strategy.}
\item{The following are the various common output options:}{}
\item{message: }{Status returned.}
\item{rowproportions: }{Vector of row proportions.}
\item{colproportions: }{Vector of column proportions.}
\item{rowclass: }{Vector of assigned row cluster to each row.}
\item{colclass: }{Vector of assigned column cluster to each column.}
\item{likelihood: }{Final pseudo log-likelihood.}
\item{rowposteriorprob: }{Final posterior probabilities for rows.}
\item{colposteriorprob: }{Final posterior probabilities for columns.}

}
}
