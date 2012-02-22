\name{cocluststrategy}
\alias{cocluststrategy}
\title{Strategy function}
\usage{
  cocluststrategy(algo = "XEMStrategy",
    initmethod = character(), stopcriteria = "Parameter",
    nbiterationsxem = 50, nbiterationsXEM = 500,
    nbinititerations = 10, initepsilon = 0.01,
    nbiterations_int = 5, epsilon_int = 0.01,
    epsilonxem = 1e-04, epsilonXEM = 1e-10, nbtry = 2,
    nbxem = 5)
}
\arguments{
  \item{algo}{The valid values for this parameter are
  "XEMStrategy" (Default) and "XCEMStrategy".This parameter
  sets the algorithm/strategy to run the model. The
  algorithms used are BEM (Block EM algorithm) for
  "XEMStrategy" and BCEM (Block classification EM
  algorithm) for "XCEMStrategy". In order to find good
  optimum, we have applied strategy to run these
  algorithms. Please refer to the documention for
  XEMStrategy for more details (You may also have to refer
  to its documentation in order to understand the
  parameters like "nbtry" and "nbxem" which plays critical
  role in obtaining good optimum).}

  \item{stopcriteria}{It specifies the stopping criteria.
  It can be based on either relative change in parameters
  vaule or relative change in log-likelihood. Valid
  criterian values are "Parameter" and "Likelihood".
  Default criteria is "Parameter".}

  \item{initmethod}{Method to initialize model paramteres.
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
  during XEM step. Deafault value is 500.}

  \item{epsilonxem}{Tolerance value used during xem step.
  Default value is 1e-4.}

  \item{epsilonXEM}{Tolerance value used during XEM step.
  Default value is 1e-10}
}
\value{
  Object of class \code{\linkS4class{strategy}}
}
\description{
  This function is used to set all the parameters for
  coclustering. It returns an object of class
  \code{\linkS4class{strategy}} which can be given as input
  to \code{\link{cocluster}} function.
}
\examples{
#Default strategy values
strategy<-cocluststrategy()
summary(strategy)
}

