#'
#' Summary function.
#' 
#' This funtion gives the summary of output from \code{cocluster}.
#' 
#' @param object output object from \code{\link{cocluster}}.
#' 
#' @name summary
#' @rdname summary-methods
#' @docType methods
#' @exportMethod summary
#' 
#' 
NULL

#' @rdname summary-methods
#' @aliases summary summary,BinaryOptions-method
#' 

setMethod(
		f="summary",
		signature = "BinaryOptions",
		definition = function(object,...) {
			if (object@successful) {
				cat("******************************************************************\n")
				cat("Model = Bernoulli Latent block model\n")
				cat("Model Parameters:")
				cat("\nClass Mean:\n")
				print(object@classmean)
				cat("\nClass Dispersion:\n")
				print(object@classdispersion)
				cat("\nRow proportions: ",object@rowproportions)
				cat("\nColumn proportions: ",object@columnproportions)
				cat("\nLikelihood: ",object@likelihood)
				cat("\n******************************************************************\n")
			} else {
				cat("Co-Clustering was not successful.\n")
			}	
		}
)


#' @rdname summary-methods
#' @aliases summary summary,ContingencyOptions-method
#'

setMethod(
		f="summary",
		signature = "ContingencyOptions",
		definition = function(object,...) {
			if (object@successful) {
				cat("******************************************************************\n")
				cat("Model Distribution = Poisson Latent block model\n")
				cat("Model Parameters:")
				cat("\nClass Gamma:\n")
				print(object@classgamma)
				cat("\nRow proportions: ",object@rowproportions)
				cat("\nColumn proportions: ",object@columnproportions)
				cat("\nLikelihood: ",object@likelihood)
				cat("\n******************************************************************\n")
			} else {
				cat("Co-Clustering was not successful.\n")
			}
		
		}
)


#' @rdname summary-methods
#' @aliases summary summary,ContinuousOptions-method
#'


setMethod(
		f="summary",
		signature = "ContinuousOptions",
		definition = function(object,...) {
			if (object@successful) {
				cat("******************************************************************\n")
				cat("Model Distribution = Gaussian Latent block model\n")
				cat("Model Parameters:")
				cat("\nClass Mean:\n")
				print(object@classmean)
				cat("\nClass Variance:\n")
				print(object@classvariance)
				cat("\nRow proportions: ",object@rowproportions)
				cat("\nColumn proportions: ",object@columnproportions)
				cat("\nLikelihood: ",object@likelihood)
				cat("\n******************************************************************\n")
			} else {
				cat("Co-Clustering was not successful.\n")
			}

		}
)
		
#' @rdname summary-methods
#' @aliases summary summary,strategy-method
#' 

setMethod(
		f="summary",
		signature = "strategy",
		definition = function(object,...) {
			cat("******************************************************************\n")
			cat("Algorithm: ",object@algo)
			cat("\nInitialization method(There is no default value): ",object@initmethod)
			cat("\nStopping Criteria: ",object@stopcriteria)
			cat("\n\nVarious Iterations")
			cat("\n******************")
			cat("\nNumber of global iterations while running initialization: ",object@nbinititerations)
			cat("\nNumber of iterations for internal E-step: ",object@nbiterations_int)
			cat("\nNumber of EM iterations used during xem: ",object@nbiterationsxem)
			cat("\nNumber of EM iterations used during XEM: ",object@nbiterationsXEM)
			cat("\nNumber of xem iterations: ",object@nbxem)
			cat("\nNumber of tries: ",object@nbtry)
			cat("\n\nVarious epsilons")
			cat("\n****************")
			cat("\nTolerance value used while initialization: ",object@initepsilon)
			cat("\nTolerance value for internal E-step: ",object@epsilon_int)
			cat("\nTolerance value used during xem: ",object@epsilonxem)
			cat("\nTolerance value used during XEM: ",object@epsilonXEM)
			
			cat("\n******************************************************************\n")
		}
)
		
