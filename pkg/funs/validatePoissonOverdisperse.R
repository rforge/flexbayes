validatePoissonOverdisperse <- function(verbose = FALSE, engine = 'OpenBUGS', n.CGR.reps=200, threshold = 0.01,
	quantilePlot.directory = NULL, plotPoorConvergence = F){
############################################################################
##
##  validatePoissonOverdisperse
##
## 	run the Cook, Gelman, and Rubin (2006) validation of a Poisson-outcome GLM.  In addition
##  to the chi-squared test to check uniformity of
##  the resulting quantiles that is specified in the paper by CGR, do an additional K-S test to
##  check the same thing.
##

	## Create the inputs to the Cook and Rubin validation
	
	# Specify the hyperparameter values.  Theoretically any values for these 
	# will lead to a valid test.  

	prior.coefficient.mean <- 						-3.5
	prior.alpha.prec <- 									2.5
	prior.beta.prec <- 										0.25
	prior.Z0 <-  													100

	hyper <- c(	prior.coefficient.mean, 
							prior.alpha.prec,
							prior.beta.prec, 
							prior.Z0  )

	# the specifications for the MCMC
	burnin <- 1000
	chainLength <- 1000
	thinning <- 20

	# the exposures for each data point
	ei <- c(1,1,1,5,5,5,5,10,10,10,10,10,10,10,20,30,90,150)
	n.data <- length (ei)

	# the inputs for the parameter generation function
	generate.param.inputs<-c(hyper, n.data)

	# the inputs for the data generation function
	generate.data.inputs<-c(n.data, ei)

	# the inputs for the function to draw from the posterior 
	analyze.data.inputs<-list(hyper, engine, chainLength, 
		thinning, burnin, n.data, ei, plotPoorConvergence)

	## run the validation
	tst <- bayes.validate(generate.param = generate.param.PoissonOverdisperse, 
		generate.param.inputs = generate.param.inputs, 
		generate.data = generate.data.PoissonOverdisperse,
		generate.data.inputs = generate.data.inputs, 
		analyze.data = analyze.data.PoissonOverdisperse, 
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


generate.param.PoissonOverdisperse <- function(inputs){
#########################################################################
##  
##  generate.param.PoissonOverdisperse
##
##  generate the parameters for the Poisson-outcome GLM from their prior distributions
##  
##  this function is used as input to the Cook and Rubin validation
##

	# get the hyperparameters
	prior.coefficient.mean <- 		inputs[1]
	alpha.prec <- 								inputs[2]
	beta.prec <- 									inputs[3]
	prior.Z0 <- 									inputs[4]
		
	# get the number of data points
	n.data <- 										inputs[5]

	# sample the precision of beta
	beta.prec <- rgamma(n = 1, shape = alpha.prec, rate = beta.prec)
	beta.standev <- 1/ sqrt(beta.prec)

	# sample beta0
	beta0 <- rnorm(n=1, mean = prior.coefficient.mean, 
			sd = beta.standev)

	# sample xi
	shrink <- runif (n=1)
	xi <- prior.Z0 * shrink / (1 - shrink)
	
	mu <- exp(beta0)
	shape <- xi
	rate <- xi / mu
	
	# sample the intensity parameters, one for each data point
	lambda <- rep(0, n.data)
	for (i in 1:n.data){
		lambda[i] <- rgamma (n=1, shape = shape, rate = rate)
	}
	
	# reparameterize to log(xi)
	logXi = log(xi)	
										
	params = c(beta0, lambda, logXi, beta.standev)
	return(params)
}

generate.data.PoissonOverdisperse <- function(params, inputs){
#########################################################################
##  
##  generate.data.PoissonOverdisperse
##
##  Given the parameters of the Poisson-outcome GLM, generate data from the model.
##  
##  This function is used as input to the Cook and Rubin validation.
##

	# get the exposures for each data point
	n.data <-  														inputs[1]
	ei <-	 																inputs[2:(n.data+1)]

	# get the parameters
	beta0 <-															params[1]
	# get the intensity parameters
	lambda <-  														params[2:(n.data+1)]
	logXi <-															params[n.data+2]
	beta.standev <- 											params[n.data+3]
			
	# sample the outcome counts
	y = rep(0,n.data)
	for (i in 1:n.data){
		poissonMean <- ei[i] * lambda[i]							
		y[i] <- rpois (n=1, lambda = poissonMean)
		}

  if (F) { #if (verbose){
  	# print the parameter values
		print("*******************************************************************")
		print("parameter values:")
		print("")
		print("beta0:")
		print(beta0)
		print("beta.standev:")
		print(beta.standev)
		print("log (xi):")
		print(logXi)
		print("intensity parameters:")
		print(lambda)	
		
		# print the data	
		print("")
		print("*****************************************************************")
		print("data:")
		print(y)
		print("")	
	}
	
	return(y)
}

analyze.data.PoissonOverdisperse <- function(y, params.true, inputs){
#########################################################################
##  
##  analyze.data.PoissonOverdisperse
##
##  Run FlexBayes on the Poisson-outcome GLM to get posterior samples.
##  
##  This function is used as input to the Cook and Rubin validation.
##

	# get the prior hyperparameters
	hyper <-														inputs[[1]]
	prior.coefficient.mean <- 					hyper[1]
	prior.alpha.prec <- 								hyper[2]
	prior.beta.prec <- 									hyper[3]
	prior.Z0 <-  												hyper[4]

	# get the chain specifications
	engine <- 													inputs[[2]]
	chainLength <- 											inputs[[3]]
	thinning <- 												inputs[[4]]
	burnin <- 													inputs[[5]]

	# get the exposures for each data point
	n.data <-	 													inputs[[6]]
	ei <- 															inputs[[7]]

	# get an indicator of whether to print
	# extra information
	plotPoorConvergence <-							inputs[[8]]

	# specify the model
	
	bugsModel <- function (){
	
	   for( i in 1 : N ) {
	      z[i] ~ dpois(poissonMean[i])
	      poissonMean[i] <- lambda[i] * e[i]
	      lambda[i] ~ dgamma(xi,invScale)
	   }
	
		# Define the precision parameter for the risk
	   xi <- (z0 * shrink) / (1 - shrink)
	   shrink ~ dunif(0,1)
		# Define the mean parameter for the risk
	   mu <- exp(beta0)
	   invScale <- xi / mu
	   beta0 ~ dnorm( beta.mean,prec)
	   prec ~ dgamma( alpha.prec, beta.prec)
	   standev <- 1 / sqrt(prec)
	   logXi <- log(xi)

 		invisible()
	}
	


	# create the data frame for fitting the model
  data <- list(z = y, e = ei, N=n.data, alpha.prec = prior.alpha.prec,
  	beta.prec = prior.beta.prec, z0 = prior.Z0, 
  	beta.mean = prior.coefficient.mean)
  
  # create the other inputs to the posteriorSamples function
	parameters.to.save <- c("beta0", "lambda", "logXi", "standev")
	initial.values <- list(shrink = 0.7, beta0 = 0, prec = 1)
	initial.values <- list(initial.values)

	posterior.samples <- posteriorSamples (data = data, model = bugsModel, 
		inits = initial.values,
		parametersToSave = parameters.to.save, nChains = 1, nThin = thinning,
		nIter = chainLength, nBurnin = burnin, DIC = F, 
		debug = F, engine = engine)
	
	# print(summary(posterior.samples))

	# check the convergence of the chain.  If not converged, print trace plots
	poorlyConvergedParams <- poorconvergenceParameters (x = posterior.samples)
	if (length(poorlyConvergedParams > 0)){
		warning(paste("Poor convergence of the following parameters: ", poorlyConvergedParams, 
			sep = ""))
		if (plotPoorConvergence)
			print(traceplot(posterior.samples[,poorlyConvergedParams]))
	}
					
	return (as.mcmc.list(posterior.samples)[[1]])
}
