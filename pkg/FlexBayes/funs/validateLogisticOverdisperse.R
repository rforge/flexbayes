validateLogisticOverdisperse <- function(verbose = FALSE, engine = 'OpenBUGS', n.CGR.reps=200, threshold = 0.01,
	quantilePlot.directory = NULL, plotPoorConvergence = F){
############################################################################
##
##  validateLogisticOverdisperse
##
## 	run the Cook, Gelman, and Rubin (2006) validation of an overdispersed 
##  logistic regression model.  
##

	if ((engine != 'WinBUGS') && (engine != 'OpenBUGS'))
		stop("engine must be one of 'WinBUGS' or 'OpenBUGS'")

	## Create the inputs to the Cook and Rubin validation
	
	# Specify the hyperparameter values.  Theoretically any values for these 
	# will lead to a valid test.  

	prior.coefficient.mean <- -10
	prior.a.prec <- 5000
	prior.b.prec <- 50
	prior.a.xi <- 5000
	prior.b.xi <- 50

	hyper <- c(	prior.coefficient.mean, 
							prior.a.prec,
							prior.b.prec, 
							prior.a.xi,
							prior.b.xi	)
			
	# the specifications for the MCMC
	burnin <- 1000
	chainLength <- 1000
	thinning <- 20

	# the number of data points and the exposures for each data point
	ni <- c(1, 5, 5, 10, 10, 30, 40, 90, 100, 130)
	n.data <- length(ni)
	
	# the inputs for the parameter generation function
	generate.param.inputs<-c(hyper, n.data)
	
	# the inputs for the data generation function
	generate.data.inputs<-c(n.data, ni, verbose)
	
	# the inputs for the function to draw from the posterior 
	analyze.data.inputs<- list(hyper, engine, chainLength, thinning, 
		burnin, n.data, ni, verbose, plotPoorConvergence)
	
	## run the validation
	tst <- bayes.validate(generate.param = generate.param.logisticOverdisperse, 
		generate.param.inputs = generate.param.inputs, 
		generate.data = generate.data.logisticOverdisperse,
		generate.data.inputs = generate.data.inputs, 
		analyze.data = analyze.data.logisticOverdisperse, 
		analyze.data.inputs = analyze.data.inputs,
		n.rep = n.CGR.reps, 
		quantilePlot.directory = quantilePlot.directory)
		
	if (verbose){
		printCGRoutput(tst, threshold)
	}
	return ((tst$adj.min.pchisq.left >= threshold) && 
		(tst$adj.min.pchisq.right >= threshold) && 
		(tst$adj.min.pnormal >= threshold))
}

generate.param.logisticOverdisperse <- function(inputs){
#########################################################################
##  
##  generate.param.logisticOverdisperse
##
##  generate the parameters for an overdispersed logistic regression 
##  model from their prior distributions
##  
##  this function is used as input to the validation methods
##

	# get the hyperparameters
	prior.coefficient.mean <- 		inputs[1]
	a.prec <- 										inputs[2]
	b.prec <- 										inputs[3]
	prior.a.xi <-  								inputs[4]
	prior.b.xi <-  								inputs[5]
	
	# get the number of data points
	n.data <- 										inputs[6]
		
	# sample the precision of beta
	prec.beta <- rgamma(n = 1, shape = a.prec, rate = b.prec)
	beta.standev <- 1/ sqrt(prec.beta)
	
	# sample beta.0
	beta.0 <- rnorm(n=1, 
		mean = prior.coefficient.mean, 
		sd = beta.standev)
	
	# sample xi
	#shrink <- runif (n=1)
	xi <- rgamma(n=1, shape=prior.a.xi, rate=prior.b.xi)
	#prior.Z0 * (shrink) / (1-shrink)
	
	# sample the intensity parameters, one for each data point
	mu <- exp(beta.0) / (1+exp(beta.0))
	shape1 <- xi * ( mu )
	shape2 <- xi * ( 1 - mu )
	lambda <- rep(0, n.data)
	for (i in (1:n.data)){
		lambda[i] <- rbeta (n=1, shape1 = shape1, shape2 = shape2)
  }
  
  params <- c(beta.0, beta.standev, lambda,  log(xi))
	return(params)
}

generate.data.logisticOverdisperse <- function(params, inputs){
#########################################################################
##  
##  generate.data.logisticOverdisperse
##
##  given the parameters of the overdispersed logistic model, 
##  generate data from the model
##  
##  this function is used as input to the validation
##
	
	# get the exposures for each data point
	n.data <-  														inputs[1]
	ni <-	 																inputs[2:(n.data+1)]
	# get an indicator of whether to print information
	verbose <- 														inputs[n.data+2]

	# get the parameters
	beta.0 <-															params[1]
	beta.standev <-  											params[2]
	log.xi <-  														params[n.data+3]
	# get the intensity parameters
	lambda <-  														params[3:(n.data+2)]
	
	
	# sample the outcome counts
	y <- rep(0,n.data)
	for (i in (1:n.data)){
		y[i] <- rbinom(n=1, size = ni[i], prob = lambda[i])
		}

  if (F) {  #if (verbose){
		# print the parameter values
		print("*******************************************************************")
		print("parameter values:")
		print("")
		print("beta.0:")
		print(beta.0)
		print("beta.standev:")
		print(beta.standev)
		print("log(xi):")
		print(log.xi)
		print("intensity parameters:")
		print(lambda)
		
		# print the data
		print("*******************************************************************")
		print("data:")
		print(y)
		print("")
	}
	return(y)
}

analyze.data.logisticOverdisperse <- function(y, params.true, inputs){
#########################################################################
##  
##  analyze.data.logisticOverdisperse
##
##  Run FlexBayes on an overdispersed logistic regression model to get 
##  posterior samples.
##  
##  This function is used as input to the validation method.
##
		
	# get the prior hyperparameters
	hyper <- 														inputs[[1]]
	prior.coefficient.mean <- 					hyper[1]
	prior.a.prec <- 										hyper[2]
	prior.b.prec <- 										hyper[3]
	prior.a.xi <-  											hyper[4]
	prior.b.xi <-  											hyper[5]
	
	# get the chain specifications
	engine <- 													inputs[[2]]
	chainLength <- 											inputs[[3]]
	thinning <- 												inputs[[4]]
	burnin <- 													inputs[[5]]
	
	# get the exposures for each data point
	n.data <-	 													inputs[[6]]
	ni <- 															inputs[[7]]
	
	# get an indicator of whether to print information
	verbose <- 						  						inputs[[8]]
	plotPoorConvergence <- 							inputs[[9]]

	# specify the model in BUGS format	
	bugsModel <- function (){
	
	   for( i in 1 : N ) {
		 	yi[i] ~ dbin (lambda[i], ni[i])
	   	lambda[i] ~ dbeta (alpha[i], beta[i])
			alpha[i] <- xi * ( mu[i])
			beta[i] <- xi * (1 - mu[i])
		 	mu[i] <- exp(meanLogitScale[i]) / (1 + exp(meanLogitScale[i]))
			meanLogitScale[i] <- beta.0	  
	  }

		# the fixed effect coefficients
	   beta.0 ~ dnorm( beta.mean, prec.beta)
	   prec.beta ~ dgamma (a.prec, b.prec)
	   beta.standev <- 1/ sqrt(prec.beta)
	   
		# define the variance parameter 
	   xi ~ dgamma(a.xi, b.xi)		
	   logXi <- log(xi)

		invisible()	
	}
	
	# create the data frame for fitting the model
	data <- list (yi = y, ni = ni, N = n.data, a.prec = prior.a.prec,
  	b.prec = prior.b.prec, a.xi = prior.a.xi, b.xi = prior.b.xi,      #z0 = prior.Z0, 
  	beta.mean = prior.coefficient.mean)
	
  # create the other inputs to the posteriorSamples function
  lambda.prior.mean <- exp(prior.coefficient.mean) / (1+ exp(prior.coefficient.mean))
	parameters.to.save <- c("beta.0", "beta.standev", "lambda", "logXi")
	initial.values = list(xi = prior.a.xi / prior.b.xi, beta.0 = prior.coefficient.mean, 
		prec.beta = prior.a.prec / prior.b.prec, lambda = rep(lambda.prior.mean, n.data))
	initial.values <- list(initial.values)

	# get the posterior samples
	posterior.samples <- posteriorSamples (data = data, model = bugsModel, 
		inits = initial.values,
		parametersToSave = parameters.to.save, nChains = 1, nThin = thinning,
		nIter = chainLength, nBurnin = burnin, DIC = F, 
		debug = F, engine = engine)
	
	# print a summary
	if (verbose){
		 # print(summary(posterior.samples))
	}

	# check the convergence of the chain.  If not converged, print trace plots and issue a warning
	poorlyConvergedParams <- poorconvergenceParameters (x = posterior.samples)
	if (length(poorlyConvergedParams > 0)){
		warning(paste("Poor convergence of the following parameters: ", poorlyConvergedParams, 
			sep = ""))
	  if (plotPoorConvergence)
			print(traceplot(posterior.samples[,poorlyConvergedParams]))
	}
					
	return (as.mcmc.list(posterior.samples)[[1]])
}
