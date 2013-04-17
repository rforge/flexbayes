validateLinearReg <- function(verbose = FALSE, engine = 'OpenBUGS', n.CGR.reps=200, threshold = 0.01, 
	quantilePlot.directory = NULL, plotPoorConvergence = F){

############################################################################
##
##  validateLinearReg
##
## 	run the Cook, Gelman, and Rubin (2006) validation of a linear regression model.  In addition
##  to the chi-squared test to check uniformity of
##  the resulting quantiles that is specified in the paper by CGR, do an additional K-S test to
##  check the same thing.
##

	## Create the inputs to the Cook and Rubin validation
	
	# Specify the hyperparameter values.  Theoretically any values for these 
	# will lead to a valid test.  

	prior.coefficient.mean <- 					-2
	prior.coefficient.prec <- 					0.3
	prior.a.prec <- 										3.5
	prior.b.prec <- 										1.5
	
	hyper <- c(	prior.coefficient.mean, 
							prior.coefficient.prec,
							prior.a.prec,
							prior.b.prec	)
	
	# the specifications for the MCMC
	burnin <- 1000
	chainLength <- 1000
	thinning <- 20

	# the number of data points and the exposures for each data point
	predictor <- c(1, 0, 5, 7.3, -2.1, 3, 12.1, 20.2, 15.7, 15.8)
	n.data <- length(predictor)
	
	# the inputs for the parameter generation function
	generate.param.inputs <- c(hyper)
	
	# the inputs for the data generation function
	generate.data.inputs <- c(n.data, predictor, verbose)
	
	# the inputs for the function to draw from the posterior 
	analyze.data.inputs <- list(hyper, engine, chainLength, 
		thinning, burnin, n.data, predictor, verbose,
		plotPoorConvergence)

	## run the validation
	tst <- bayes.validate(generate.param = generate.param.linearReg, 
		generate.param.inputs = generate.param.inputs, 
		generate.data = generate.data.linearReg,
		generate.data.inputs = generate.data.inputs, 
		analyze.data = analyze.data.linearReg, 
		analyze.data.inputs = analyze.data.inputs,
		n.rep = n.CGR.reps, quantilePlot.directory = quantilePlot.directory)
  
	if (verbose){
		printCGRoutput(tst, threshold)
	}
	return ((tst$adj.min.pchisq.left >= threshold) && 
		(tst$adj.min.pchisq.right >= threshold) && 
		(tst$adj.min.pnormal >= threshold))
}

generate.param.linearReg <- function(inputs){
#########################################################################
##  
##  generate.param.linearReg
##
##  generate the parameters for the linear regression model from their prior distributions
##  
##  this function is used as input to the Cook and Rubin validation
##

	# get the hyperparameters
	prior.coefficient.mean <- 						inputs[1]
	prior.coefficient.prec <- 						inputs[2]
	a.prec <- 														inputs[3]
	b.prec <- 														inputs[4]
		
	# sample the beta parameters
	prior.coef.sd <- 1/ (prior.coefficient.prec)^0.5
	beta1 <- rnorm(n=1, 
		mean = prior.coefficient.mean, 
		sd = prior.coef.sd)
	beta2 <- rnorm(n=1, 
		mean = prior.coefficient.mean, 
		sd = prior.coef.sd)
		
	# sample the precision
	tau.c <- rgamma(n = 1, shape = a.prec, rate = b.prec)
	sd <- 1/sqrt( tau.c )
		  
  params <- c(beta1, beta2, sd)
	return(params)
}

generate.data.linearReg <- function(params, inputs){
#########################################################################
##  
##  generate.data.linearReg
##
##  given the parameters of the linear regression model, generate data from the model
##  
##  this function is used as input to the Cook and Rubin validation
##
	
	# get the number of data points
	n.data <-  														inputs[1]
	# get the predictor for each data point
	x <-	 																inputs[2:(n.data+1)]
	# get an indicator of whether to print information
	verbose <- 														inputs[n.data+2]
	# get the parameters
	beta1 <-  														params[1]
	beta2 <-															params[2]
	sd <-  														params[3]
	
	# sample the outcome
	y <- rep(0,n.data)
	for (i in (1:n.data)){
		mu <- beta1 + beta2*x[i]
		y[i] <- rnorm(n=1, mean = mu, sd = sd)
		}
	
	if (F) { #if (verbose){
		# print the parameter values
		print("*******************************************************************")
		print("parameter values:")
		print("")
		print("beta1:")
		print(beta1)
		print("beta2:")
		print(beta2)
		print("sd")
		print(sd)
	
		# print the data
		print("*******************************************************************")
		print("data:")
		print(y)
		print("")
	}
	
	return(y)
}

analyze.data.linearReg <- function(y, params.true, inputs){
#########################################################################
##  
##  analyze.data.linearReg
##
##  run Splus Bayes on the linear regression model to get posterior samples
##  
##  this function is used as input to the Cook and Rubin validation
##
		
	# get the prior hyperparameters
	hyper <-														inputs[[1]]
	prior.coefficient.mean <- 					hyper[1]
	prior.coefficient.prec <- 					hyper[2]
	prior.a.prec <- 										hyper[3]
	prior.b.prec <- 										hyper[4]
	
	# get the chain specifications
	engine <-     											inputs[[2]]
	chainLength <- 											inputs[[3]]
	thinning <- 												inputs[[4]]
	burnin <- 													inputs[[5]]
	
	# get the exposures for each data point
	n.data <-	 													inputs[[6]]
	x <- 																inputs[[7]] 
	verbose <- 													inputs[[8]] 
	plotPoorConvergence <- 							inputs[[9]]
	
	if ( is.element( engine, c("WinBUGS", "OpenBUGS") ) ){
		# specify the model in BUGS format
		bugsModel <- function (){
			for( i in 1 : N ) {
	
				outcome[i] ~ dnorm(mu[i],tau.c)
				mu[i] <- beta1 + beta2 * predictor[i]
			}
			beta1 ~ dnorm(beta.mean, beta.prec)
			beta2 ~ dnorm(beta.mean, beta.prec)
	
			tau.c ~ dgamma(a.prec, b.prec)
			sigma <- 1 / sqrt(tau.c)
	
			invisible()		
		}
		# specify the model in BUGS format	
		
		# create the data frame for fitting the model
		data <- list (outcome = y, predictor = x, N = n.data, a.prec = prior.a.prec,
	  	b.prec = prior.b.prec, 
	  	beta.mean = prior.coefficient.mean, beta.prec = prior.coefficient.prec)
		
	  # create the other inputs to the posteriorSamples function
		parameters.to.save <- c("beta1", "beta2", "sigma")
		initial.values <- list(beta1 = prior.coefficient.mean, 
			beta2 = prior.coefficient.mean, 
			tau.c = prior.a.prec / prior.b.prec)
		initial.values <- list(initial.values)
	
		# get the posterior samples
		posterior.samples <- posteriorSamples (data = data, model = bugsModel, 
			parametersToSave = parameters.to.save, nChains = 1, nThin = thinning,
			nIter = chainLength, nBurnin = burnin, DIC = F, inits = initial.values,
			debug = F, engine = engine)
	} else {
	
		linear.sampler <- blm.sampler( nThin = thinning, nSamples = chainLength,
			nBurnin = burnin )
		
		beta.prior <- bayes.normal( rep(prior.coefficient.mean, 2), 
			diag( rep ( 1/prior.coefficient.prec, 2 ) ) )
			
		sigma.prior <- bayes.invChisq( df=(2*prior.a.prec), 
			sigma0.sq=(prior.b.prec/prior.a.prec) )
			
		linear.prior <- blm.prior ( fixed.coef = beta.prior, sigma2 = 
			 sigma.prior )
  		
		posterior.samples <- blm ( data = data.frame( y=y, x=x ), 
			prior = linear.prior, sampler = linear.sampler, 
			formula = y ~ x, likelihood = blm.likelihood () )
	}
	if (verbose)	{
		# print a summary
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
