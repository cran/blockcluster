#'
#' @include optionclasses.R
#' @include Strategy.R
#' @include summary.R
#' @include plot.R
#' @include onload.R
#' @include Rcoclust.R
#' 
NULL


#' Co-Clustering function.
#' 
#' This function performs Co-Clustering (simultanuous clustering of rows and columns ) for Binary, Contingency
#' and Continuous datasets using latent block models.  
#' 
#' @param data Input data as matrix (or list containing data matrix, numeric vector for row effects and numeric 
#'        vector column effects in case of contingency data with known row and column effects.)
#' @param datatype This is the type of data which can be "binary" , "contingency" or "continuous".
#' @param model This is the name of model. The following models exists for various kinds of datasets:
#' \tabular{rlll}{
#'     Model  \tab Datatype \tab Proportions \tab Dispersion/Variance \cr
#'     pik_rhol_epsilonkl(Default) \tab binary \tab unequal \tab unequal \cr
#'     pik_rhol_epsilon \tab binary \tab unequal \tab equal \cr
#'     pi_rho_epsilonkl \tab binary \tab equal \tab unequal \cr
#'     pi_rho_epsilon \tab binary \tab equal \tab equal \cr
#'     pik_rhol_sigma2kl \tab continuous \tab unequal \tab unequal \cr
#'     pik_rhol_sigma \tab continuous \tab unequal \tab equal \cr
#'     pi_rho_sigma2kl \tab continuous \tab equal \tab unequal \cr
#'     pi_rho_sigma2 \tab continuous \tab equal \tab equal \cr
#'     pik_rhol_unknown(default) \tab contingency \tab unequal \tab N.A \cr
#'     pi_rho_unknown \tab contingency \tab equal \tab N.A \cr
#'     pik_rhol_known \tab contingency \tab unequal \tab N.A \cr
#'     pi_rho_known \tab contingency \tab equal \tab N.A \cr
#' }
#' 
#' @param nbcocluster Interger vector specifying the number of row and column clusters respectively.
#' @param strategy Object of class \code{\linkS4class{strategy}}.
#' @return Return an object of \code{\linkS4class{BinaryOptions}} or \code{\linkS4class{ContingencyOptions}}
#' or \code{\linkS4class{ContinuousOptions}} depending of whether the datatype is Bianry, Contingency or Continuous
#' respectively.
#' 
#' @export
#' 
#' @exportPattern "^[[:alpha:]]+"
#' @useDynLib blockcluster
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
#' #Plot the original and co-clustered data 
#' plot(out)
#' 
#' # A little advanced example with simulated gaussian data
#' data(gaussiandata)
#' #set strategy , see documentation for "cocluststrategy" function for more details.
#' newstrategy<-cocluststrategy(nbxem=5,nbtry=2,algo="XEMStrategy")
#' # calling cocluster function with newstrategy and default model
#' out<-cocluster(gaussiandata,datatype="continuous",nbcocluster=c(2,3),strategy=newstrategy)
#' 

cocluster<-function(data, datatype, model, nbcocluster, strategy = cocluststrategy()) 
{
	#Check for data
	if(missing(data)){
		stop("Data is missing.")
	}
	else
	{
		if(!is.list(data)){
			if(!is.matrix(data))
				stop("Data should be matrix.")
		}else
		{
			if(!is.matrix(data[[1]]))
				stop("Data should be matrix.")
			if(!is.numeric(data[[2]])||!is.numeric(data[[3]]))
				stop("Row/Column effects should be numeric.")
			if(length(data[[2]])!=dim(data[[1]])[1]||length(data[[3]])!=dim(data[[1]])[2])
				stop("Dimention mismatch in Row/column effects  and Data.")
		}
	}
	
	#check for number of coclusters
	if(missing(nbcocluster))
	{
		stop("Mention number of CoClusters.")
	}
	else
	{
		if(!is.list(data))
		dimention = dim(data)
		else
			dimention = dim(data[[1]])
		
		if(dimention[1]<nbcocluster[1])
			stop("Number of Row cluters exceeds numbers of rows.")
		
		if(dimention[2]<nbcocluster[2])
			stop("Number of Column cluters exceeds numbers of columns.")
	}
	
	#check for Algorithm name
	
	if(strategy@algo!="XEMStrategy" && strategy@algo!="XCEMStrategy")
		stop("Incorrect Algorithm, Valide algorithms are: XEMStrategy, XCEMStrategy") else{}
	
	#check for stopping criteria
	
	if(strategy@stopcriteria!="Parameter" && strategy@stopcriteria!="Likelihood")
		stop("Incorrect stopping criteria, Valid stopping criterians are: Parameter, Likelihood")else{}
	
	#check for datatype and models
	if (missing(datatype)) {
		stop("Mention datatype.")
	} 
	else{
		if(datatype == "binary"){
			if(missing(model)){
				model = "pik_rhol_epsilonkl"
			}
			else 
			{
				if(model!="pik_rhol_epsilonkl" && model!="pik_rhol_epsilon" && 
						model!="pi_rho_epsilonkl" && model!="pi_rho_epsilon")
				{
					stop("Incorrect Model, Valid Binary models are:pik_rhol_epsilonkl, pik_rhol_epsilon
									pi_rho_epsilonkl, pi_rho_epsilon")
				}
			}
			
			if(length(strategy@initmethod)==0){
				strategy@initmethod = "CEMInit"
			}
			else
			{
				if(strategy@initmethod!="CEMInit")
					stop("Incorrect initialization method, valid method(s) are: CEMInit")
			}
			
			inpobj<-new("BinaryOptions",data = data,datatype = datatype, model = model,nbcocluster = nbcocluster, strategy = strategy)
			
		}
		
		else if(datatype == "continuous"){
			if(missing(model)){
				model = "pik_rhol_sigma2kl"
			}
			else if(model!="pik_rhol_sigma2kl" && model!="pik_rhol_sigma2" && 
					model!="pi_rho_sigma2kl" && model!="pi_rho_sigma2"){
				stop("Incorrect Model, Valid Binary models are: pik_rhol_sigma2kl, pik_rhol_sigma2, pi_rho_sigma2kl, pi_rho_sigma2")
			}
			
			if(length(strategy@initmethod)==0){
				strategy@initmethod = "CEMInit"
			}
			else
			{
				if(strategy@initmethod!="CEMInit")
					stop("Incorrect initialization method, valid method(s) are: CEMInit")
			}
			
			inpobj<-new("ContinuousOptions",data = data,datatype = datatype, model = model,nbcocluster = nbcocluster,strategy = strategy)
			
		}
		
		else if(datatype == "contingency"){
			if(missing(model)&& !is.list(data)){
				model = "pik_rhol_unknown"
			}
			else if(missing(model) && is.list(data))
			{
				model = "pik_rhol_known"
			}
			else if(model!="pik_rhol_unknown" && model!="pik_rhol_known" && 
					model!="pi_rho_unknown" && model!="pi_rho_known"){
				stop("Incorrect Model, Valid Binary models are:pik_rhol_unknown, pik_rhol_known, pi_rho_unknown, pi_rho_known")
		}
		
		if((model=="pi_rho_known"||model=="pik_rhol_known")&& (length(data)!=3))
		{
			stop("Missing Row/Column effects.") 
		}
		
		if(length(strategy@initmethod)==0){
			if((model=="pi_rho_known"||model=="pik_rhol_known"))
			{strategy@initmethod = "RandomInit"}
			else{strategy@initmethod = "CEMInit"}
		}
		else
		{
			if(strategy@initmethod!="RandomInit"&&(model=="pi_rho_known"||model=="pik_rhol_known"))
			{stop("Incorrect initialization method, valid method(s) are: RandomInit")}
			else if(strategy@initmethod!="CEMInit"&&(model=="pi_rho_unknown"||model=="pik_rhol_unknown"))
				stop("Incorrect initialization method, valid method(s) are: CEMInit")
		}
		if(!is.list(data))
		inpobj<-new("ContingencyOptions",data = data,datatype = datatype, model = model, nbcocluster = nbcocluster, strategy = strategy)
		else
			inpobj<-new("ContingencyOptions",data = data[[1]],datatype = datatype, model = model, 
					nbcocluster = nbcocluster, strategy = strategy,datamui=data[[2]],datanuj=data[[3]])
		
	}
		else
		{stop("Invalid datatype, Valid types are: binary , contingency , continuous")}
	}
	

	
   .Call("CoClustmain",inpobj,PACKAGE = "blockcluster")
	
   print(inpobj@message)
   
   return(inpobj)
}

