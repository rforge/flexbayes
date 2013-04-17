validatePoissonGLM <- function(verbose = FALSE, engine = 'OpenBUGS', n.CGR.reps=200, threshold = 0.01,
	quantilePlot.directory = NULL, plotPoorConvergence = F){
############################################################################
##
##  validatePoissonGLM
##
## 	run the Cook, Gelman, and Rubin (2006) validation of a Poisson-outcome 
##  GLM. 
##

	if ( !is.element( engine, c( 'WinBUGS', 'OpenBUGS', 'bpm' ) ) )
		stop("engine must be one of 'WinBUGS', 'OpenBUGS', or 'bpm'")

	## Create the inputs to the Cook and Rubin validation
	
	# Specify the hyperparameter values.  Theoretically any values for these 
	# will lead to a valid test.  

	prior.coef.mean <- 						-0.5
	prior.coef.var <- 						0.7
	
	hyper <- c(	prior.coef.mean, 
		prior.coef.var )

	# the specifications for the MCMC
	burnin <- 1000
	chainLength <- 1000
	thinning <- 20

	# the exposures for each data point
	ei <- c(1,1,1,5,5,5,5,10,10,10,10,10,10,10,20,30,90,150)
	n.data <- length (ei)
	predictor <- c( rep (0.2, 8), rep(-0.3, 5), rep( 0.3, 5 ) )

	# the inputs for the parameter generation function
	generate.param.inputs <- hyper

	# the inputs for the data generation function
	generate.data.inputs <- list(n.data, ei, predictor)

	# the inputs for the function to draw from the posterior 
	analyze.data.inputs<-list(hyper, engine, chainLength, 
		thinning, burnin, n.data, ei, predictor, 
		plotPoorConvergence)

	## run the validation
	tst <- bayes.validate(generate.param = generate.param.PoissonGLM, 
		generate.param.inputs = generate.param.inputs, 
		generate.data = generate.data.PoissonGLM,
		generate.data.inputs = generate.data.inputs, 
		analyze.data = analyze.data.PoissonGLM, 
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

generate.param.PoissonGLM <- function(inputs){
#########################################################################
##  
##  generate.param.PoissonGLM
##
##  generate the parameters for the Poisson-outcome GLM from their prior 
##  distributions.
##
##  this function is used as input to the Cook and Rubin validation
##

	# get the hyperparameters
	prior.coef.mean <- inputs[1]
	prior.coef.var <- inputs[2]
		
	# sample beta0 and beta1
	sd <- sqrt(prior.coef.var)
	beta0 <- rnorm(n=1, mean = prior.coef.mean, 
			sd = sd)
	beta1 <- rnorm(n=1, mean = prior.coef.mean, 
			sd = sd)

	params = c(beta0, beta1)
	return(params)
}

generate.data.PoissonGLM <- function(params, inputs){
#########################################################################
##  
##  generate.data.PoissonGLM
##
##  Given the parameters of the Poisson-outcome GLM, generate data from 
##  the model.
##  This function is used as input to the Cook and Rubin validation.
##

	# get the exposures for each data point
	n.data <- inputs[[1]]
	ei <-	inputs[[2]]
	predictor <- inputs[[3]]

	# sample the outcome counts
	y <- rep(0,n.data)
	beta <- params
	for (i in 1:n.data){
		lambda <- exp ( beta[1] + beta[2] * predictor[i] )
		poissonMean <- ei[i] * lambda
		y[i] <- rpois (n=1, lambda = poissonMean)
		}

  if (F) { #if (verbose){
  	# print the parameter values
		print("*******************************************************************")
		print("parameter values:")
		print("")
		print("beta0:")
		print(beta[1])
		print("beta1:")
		print(beta[2])
		
		# print the data	
		print("")
		print("*****************************************************************")
		print("data:")
		print(y)
		print("")	
	}
	
	return(y)
}

analyze.data.PoissonGLM <- function(y, params.true, inputs){
#########################################################################
##  
##  analyze.data.PoissonGLM
##
##  Run FlexBayes on the Poisson-outcome GLM to get posterior samples.
##  
##  This function is used as input to the Cook and Rubin validation.
##

	# get the prior hyperparameters
	hyper <- inputs[[1]]
	prior.coef.mean <- hyper[1]
	prior.coef.var <- hyper[2]
	
	# get the chain specifications
	engine <- inputs[[2]]
	chainLength <- inputs[[3]]
	thinning <- inputs[[4]]
	burnin <- inputs[[5]]

	# get the exposures for each data point
	n.data <-	inputs[[6]]
	ei <- inputs[[7]]
  predictor <- inputs[[8]]
	
	# get an indicator of whether to print
	# extra information
	plotPoorConvergence <-							inputs[[9]]

	if (is.element(engine, c("OpenBUGS", "WinBUGS"))){
	
		# specify the model	
		bugsModel <- function (){
		
		   for( i in 1 : N ) {
		      y[i] ~ dpois(poissonMean[i])
		      poissonMean[i] <- lambda[i] * e[i]
		      lambda[i] <- exp ( beta0 + beta1 * predictor[i] )
		   }
		
		   beta0 ~ dnorm( beta.mean, beta.prec)
		   beta1 ~ dnorm( beta.mean, beta.prec)
		   
	 		invisible()
		}
		
		# create the data frame for fitting the model
	  data <- list(y = y, e = ei, predictor = predictor,
	  	N=n.data, beta.mean = prior.coef.mean,
	  	beta.prec = 1/(prior.coef.var) )
	  
	  # create the other inputs to the posteriorSamples function
		parameters.to.save <- c("beta0", "beta1")
		initial.values <- list(beta0 = prior.coef.mean, beta1 = prior.coef.mean)
		initial.values <- list(initial.values)
	
		posterior.samples <- posteriorSamples (data = data, model = bugsModel, 
			inits = initial.values,
			parametersToSave = parameters.to.save, nChains = 1, nThin = thinning,
			nIter = chainLength, nBurnin = burnin, DIC = F, 
			engine = engine)
	} else {
		
		beta.prior <- bayes.normal( rep( prior.coef.mean, 2 ), 
				diag( rep( prior.coef.var, 2 ) ) )
		poisson.prior <- bpm.prior( fixed.coef = beta.prior )
		
		poisson.sampler <- bpm.sampler( nSamples = chainLength, nBurnin = burnin,
			nThin = thinning, nChains = 1 )
			
		data <- data.frame( y=y, e=ei, predictor=predictor )
			
		posterior.samples <- bpm( data = data,
			exposure.formula = ~ e, fixed.formula = y ~ predictor, 
			prior = poisson.prior, sampler = poisson.sampler )
		
	}
	
	# print(summary(posterior.samples))

	# check the convergence of the chain.  If not converged, print trace plots
	poorlyConvergedParams <- poorconvergenceParameters (x = posterior.samples)
	if (length(poorlyConvergedParams > 0)){
		warning(paste("Poor convergence of the following parameters: ", 
		  poorlyConvergedParams, sep = ""))
		if (plotPoorConvergence)
			print(traceplot(posterior.samples[,poorlyConvergedParams]))
	}
					
	return (as.mcmc.list(posterior.samples)[[1]])
}
