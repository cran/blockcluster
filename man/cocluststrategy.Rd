\name{cocluststrategy}
\alias{cocluststrategy}
\title{Strategy function}
\usage{
  cocluststrategy(algo = "BEM", initmethod = character(),
    stopcriteria = "Parameter", nbiterationsxem = 50,
    nbiterationsXEM = 500, nbinititerations = 10,
    initepsilon = 0.01, nbiterations_int = 5,
    epsilon_int = 0.01, epsilonxem = 1e-04,
    epsilonXEM = 1e-10, nbtry = 2, nbxem = 5,
    bayesianform = FALSE, hyperparam = numeric(0))
}
\arguments{
  \item{algo}{The valid values for this parameter are "BEM"
  (Default), "BCEM" and "BSEM".}

  \item{stopcriteria}{It specifies the stopping criteria.
  It can be based on either relative change in parameters
  (preffered due to computation reasons) value or relative
  change in pseudo log-likelihood. Valid criterion values
  are "Parameter" and "Likelihood". Default criteria is
  "Parameter".}

  \item{initmethod}{Method to initialize model parameters.
  The valid values are "CEMInit", "FuzzyCEMInit" and
  "RandomInit". For now only one kind of initialization
  exist for every model currently available in the package.
  Hence default value for initialization is set according
  to the model.}

  \item{nbinititerations}{Number of Global iterations used
  in initialization step. Default value is 10.}

  \item{initepsilon}{Tolerance value used while
  initialization. Default value is 1e-2.}

  \item{nbiterations_int}{Number of iterations for internal
  E step. Default value is 5.}

  \item{epsilon_int}{Tolerance value for relative change in
  Parameter/likelihood for internal E-step. Default value
  is 1e-2.}

  \item{nbtry}{Number of tries (XEM steps). Default value
  is 2.}

  \item{nbxem}{Number of xem steps. Default value is 5.}

  \item{nbiterationsxem}{Number of EM iterations used
  during xem step. Default value is 50.}

  \item{nbiterationsXEM}{Number of EM iterations used
  during XEM step. Default value is 500.}

  \item{epsilonxem}{Tolerance value used during xem step.
  Default value is 1e-4.}

  \item{epsilonXEM}{Tolerance value used during XEM step.
  Default value is 1e-10}

  \item{bayesianform}{Boolean parameter to indicate whether
  to run algorithms in bayesian settings or not.Default
  value is false.}

  \item{hyperparam}{Hyper-parameters ("a" and "b") in case
  of Bayesian settings.}
}
\value{
  Object of class \code{\linkS4class{strategy}}
}
\description{
  This function is used to set all the parameters for
  Co-clustering. It returns an object of class
  \code{\linkS4class{strategy}} which can be given as input
  to \code{\link{cocluster}} function.
}
\examples{
#Default strategy values

strategy<-cocluststrategy()
summary(strategy)
}

