validateLogisticTprior <- function(verbose = FALSE, engine = 'OpenBUGS', 
  n.CGR.reps=200, threshold = 0.01,
	quantilePlot.directory = NULL, plotPoorConvergence = F){
############################################################################
##
##  validateLogisticTprior
##
## 	run the Cook, Gelman, and Rubin (2006) validation of a logistic 
##  regression model with t prior distribution on the coefficients.  
##

	if ( !is.element( engine, c( 'WinBUGS', 'OpenBUGS', 'bbm' ) ) )
		stop("engine must be one of 'WinBUGS', 'OpenBUGS', or 'bbm'")

	## Create the inputs to the Cook and Rubin validation
	
	# Specify the hyperparameter values.  Theoretically any values for these 
	# will lead to a valid test.  

	prior.coef.mean <- c( -1.5, 0.3)
	prior.coef.cov <- matrix( c( 3.2, 0.01, 0.01, 2.5 ), nrow = 2 )
	prior.coef.df <- 4.2

	hyper <- list( prior.coef.mean, 
							prior.coef.cov, prior.coef.df )
			
	# the specifications for the MCMC
	burnin <- 1000
	chainLength <- 1000
	thinning <- 20
	
	# the number of data points
	n.data <- 10
	trials <- c( rep(8, 5), rep(20, 5) )
	predictor <- c( rep( 0.2, 3 ), rep( 0.7, 2), rep ( -1.1, 5 ) )

	# the inputs for the parameter generation function
	generate.param.inputs <- hyper
	
	# the inputs for the data generation function
	generate.data.inputs <- list(n.data, trials, predictor, verbose)
	
	# the inputs for the function to draw from the posterior 
	analyze.data.inputs<- list(hyper, engine, chainLength, thinning, 
		burnin, n.data, trials, predictor, verbose, plotPoorConvergence)
	
	## run the validation
	tst <- bayes.validate(generate.param = generate.param.logisticTprior, 
		generate.param.inputs = generate.param.inputs, 
		generate.data = generate.data.logisticTprior,
		generate.data.inputs = generate.data.inputs, 
		analyze.data = analyze.data.logisticTprior, 
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

rmvt <- function(mean, cov, df ){
#########################################################################
##  
##  rmvt
##
##  generate from the multivariate t distribution 
##  

	# get the Cholesky decomposition of the covariance matrix
	A <- chol( cov )

	x <- rchisq( n=1, df=df )
	z <- rnorm( n = length(mean) )
	z <- matrix( z, ncol = 1 )
	tRV <- A %*% z 

	tRV <- tRV * sqrt( df / x )
	tRV <- tRV + mean
	
	return( as.vector( tRV ) )
}

generate.param.logisticTprior <- function(inputs){
#########################################################################
##  
##  generate.param.logisticTprior
##
##  generate the parameters for a logistic regression model from their 
##  prior distributions
##  
##  this function is used as input to the Cook and Rubin validation
##

	# get the hyperparameters
	prior.coef.mean <- inputs[[1]]
	prior.coef.cov <- inputs[[2]]
	prior.coef.df <- inputs[[3]]

	# sample beta
	beta <- rmvt( mean = prior.coef.mean, cov = prior.coef.cov, 
	  df = prior.coef.df )

	return(beta)
}

generate.data.logisticTprior <- function(params, inputs){
#########################################################################
##  
##  generate.data.logisticTprior
##
##  given the parameters of a logistic regression model, generate data 
##  from the model
##  
##  this function is used as input to the Cook and Rubin validation
##
	
	# get the exposures for each data point
	n.data <-  														inputs[[1]]
	trials <-															inputs[[2]]
	predictor <-													inputs[[3]]
	# get an indicator of whether to print information
	verbose <- 														inputs[[4]]

	# get the parameters
	beta <-															params
	
	# sample the outcome counts
	y <- rep(0,n.data)
	for (i in (1:n.data)){
		mean <- beta[1] + beta[2] * predictor[i]
		prob <- exp(mean)/(1+exp(mean))
		y[i] <- rbinom(n=1, size = trials[i], prob = prob)
		}

  if (F) { #if (verbose){
		# print the parameter values
		print("*******************************************************************")
		print("parameter values:")
		print("")
		print("beta")
		print(beta)
		
		# print the data
		print("*******************************************************************")
		print("data:")
		print(y)
		print("")
	}
	return(y)
}

analyze.data.logisticTprior <- function(y, params.true, inputs){
#########################################################################
##  
##  analyze.data.logisticTprior
##
##  Run FlexBayes on a logistic regression model to get posterior samples.
##  
##  This function is used as input to the Cook and Rubin validation.
##
		
	# get the prior hyperparameters
	hyper <- 														inputs[[1]]
	prior.coef.mean <- hyper[[1]]
	prior.coef.cov <- hyper[[2]]
	prior.coef.df <- hyper[[3]]
	
	# get the chain specifications
	engine <- 													inputs[[2]]
	chainLength <- 											inputs[[3]]
	thinning <- 												inputs[[4]]
	burnin <- 													inputs[[5]]
	
	# get the exposures for each data point
	n.data <-	inputs[[6]]
	trials <- inputs[[7]]
	predictor <- inputs[[8]]
	
	# get an indicator of whether to print information
	verbose <- 						  						inputs[[9]]
	plotPoorConvergence <- 							inputs[[10]]

	if ( is.element( engine, c("OpenBUGS", "WinBUGS") ) ){
		# specify the model in BUGS format	
		bugsModel <- function (){
		
		   for( i in 1 : N ) {
			 	y[i] ~ dbin (lambda[i], trials[i])
		   	logit(lambda[i]) <- mean[i]
		   	mean[i] <- beta[1] + beta[2] * predictor[i]
		  }
	
			# the fixed effect coefficients
		  beta[1:2] ~ dmt( beta.mean[], beta.prec[,], beta.df )
			
			invisible()	
		}
		
		# create the data frame for fitting the model
		data <- list (y = y, trials = trials, 
			predictor = predictor, N = n.data, 
	  	beta.mean = prior.coef.mean,
	  	beta.prec = solve( prior.coef.cov ),
	  	beta.df = prior.coef.df )
		
	  # create the other inputs to the posteriorSamples function
		parameters.to.save <- c("beta")
		initial.values = list( beta = prior.coef.mean )
		initial.values <- list(initial.values)
	
		# get the posterior samples
		posterior.samples <- posteriorSamples (data = data, model = bugsModel, 
			inits = initial.values,
			parametersToSave = parameters.to.save, nChains = 1, nThin = thinning,
			nIter = chainLength, nBurnin = burnin, DIC = F, 
			debug = F, engine = engine)
	} else {
	
		beta.prior <- bayes.t( mean.vector = prior.coef.mean, 
				covmat = prior.coef.cov, df = prior.coef.df )
		logistic.prior <- bpm.prior( fixed.coef = beta.prior )
		
		logistic.sampler <- bpm.sampler( nSamples = chainLength, nBurnin = burnin,
			nThin = thinning, nChains = 1 )
			
		data <- data.frame( y=y, trials=trials, predictor=predictor )
			
		posterior.samples <- bbm( data = data,
			trials.formula = ~ trials, fixed.formula = y ~ predictor, 
			prior = logistic.prior, sampler = logistic.sampler,
			engine = "WinBUGS" )
		
	}
	# print a summary
	if (verbose){
		 # print(summary(posterior.samples))
	}

	# check the convergence of the chain.  If not converged, print trace plots 
	# and issue a warning
	poorlyConvergedParams <- poorconvergenceParameters (x = posterior.samples)
	if (length(poorlyConvergedParams > 0)){
		warning(paste("Poor convergence of the following parameters: ", 
			poorlyConvergedParams, sep = ""))
	  if (plotPoorConvergence)
			print(traceplot(posterior.samples[,poorlyConvergedParams]))
	}
					
	return (as.mcmc.list(posterior.samples)[[1]])
}
