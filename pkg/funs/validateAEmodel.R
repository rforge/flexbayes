validateAEmodel <- function(verbose = F, engine = 'OpenBUGS', n.CGR.reps=200, threshold = 0.01, 
	quantilePlot.directory = NULL, plotPoorConvergence = F){

############################################################################
##
##  validateAEmodel
##
## 	run the Cook, Gelman, and Rubin (2006) validation of the adverse event
##  model given in Berry and Berry (2004).  In addition to the chi-squared 
##  test to check uniformity of the resulting quantiles that is specified 
##  in the paper by CGR, do an additional K-S test to check the same thing.
##

	## Create the inputs to the Cook and Rubin validation
	
	# Specify B, Nobs, Nc, and Nt
	# B: # of body systems
	# Nobs: # of AE types in each body system
	# Nc: # of individuals in control group
	# Nt: # of individuals in treatment group
	B <- 3
	Nobs <- rep(5, B)
	Nc <- 50
	Nt <- 50

	# Specify the hyperparameter values.  Theoretically any values for these 
	# will lead to a valid test.  

	hyper <- list(	
		lambda.alpha = 1, lambda.beta = 1, 
		alpha.sigma.theta = 3, beta.sigma.theta = 1,
		alpha.sigma.gamma = 3, beta.sigma.gamma = 1, 
		alpha.tau.theta = 3, beta.tau.theta = 1,
		alpha.tau.gamma = 3, beta.tau.gamma = 1, 
		mu.gamma.00 = 0, tau2inv.gamma.00 = 1,
		mu.theta.00 = 0, tau2inv.theta.00 = 1 )
	
	# the specifications for the MCMC
	burnin <- 10000
	chainLength <- 1000
	thinning <- 10
	
	# the inputs for the parameter generation function
	generate.param.inputs <- list(B, Nobs, Nc, Nt, hyper, verbose)
	
	# the inputs for the data generation function
	generate.data.inputs <- list(B, Nobs, Nc, Nt, verbose)
	
	# the inputs for the function to draw from the posterior 
	analyze.data.inputs <- list(hyper, engine, chainLength, 
		thinning, burnin, verbose,
		plotPoorConvergence)

	## run the validation
	tst <- bayes.validate(generate.param = generate.param.AEmodel, 
		generate.param.inputs = generate.param.inputs, 
		generate.data = generate.data.AEmodel,
		generate.data.inputs = generate.data.inputs, 
		analyze.data = analyze.data.AEmodel, 
		analyze.data.inputs = analyze.data.inputs,
		n.rep = n.CGR.reps, quantilePlot.directory = quantilePlot.directory)
  
	if (verbose){
		printCGRoutput(tst, threshold)
	}
	return ((tst$adj.min.pchisq.left >= threshold) && 
		(tst$adj.min.pchisq.right >= threshold) && 
		(tst$adj.min.pnormal >= threshold))
}


generate.param.AEmodel <- function(inputs){
#########################################################################
##  
##  generate.param.AEmodel
##
##  generate the parameters from their prior distributions
##  
##  this function is used as input to the Cook and Rubin validation
##

	# get B, Nobs, Nc, and Nt
	B <- inputs[[1]]
	Nobs <- inputs[[2]]
	Nc <- inputs[[3]]
	Nt <- inputs[[4]]

	# get the hyperparameters
	hyper <- inputs[[5]]
	lambda.alpha <- hyper$lambda.alpha
	lambda.beta <- hyper$lambda.beta
	alpha.sigma.theta <- hyper$alpha.sigma.theta
	beta.sigma.theta <- hyper$beta.sigma.theta
	alpha.sigma.gamma <- hyper$alpha.sigma.gamma
	beta.sigma.gamma <- hyper$beta.sigma.gamma
	alpha.tau.theta <- hyper$alpha.tau.theta
	beta.tau.theta <- hyper$beta.tau.theta
	alpha.tau.gamma <- hyper$alpha.tau.gamma
	beta.tau.gamma <- hyper$beta.tau.gamma
		
	mu.gamma.00 <- hyper$mu.gamma.00
	tau2inv.gamma.00 <- hyper$tau2inv.gamma.00
	mu.theta.00 <- hyper$mu.theta.00
	tau2inv.theta.00 <- hyper$tau2inv.theta.00

	# Get an indicator of whether to print information
	verbose <- inputs[[6]]

	# Sample alpha.pi and beta.pi
	alphaExp <- rexp(n=1, rate=lambda.alpha)
	alpha.pi <- alphaExp + 1
	betaExp <- rexp(n=1, rate=lambda.beta)
	beta.pi <- betaExp + 1

	# Sample mu.gamma.0 and mu.theta.0
	mu.gamma.0 <- rnorm(n=1, mean = mu.gamma.00, 
		sd = 1/sqrt(tau2inv.gamma.00))
	mu.theta.0 <- rnorm(n=1, mean = mu.theta.00, 
		sd = 1/sqrt(tau2inv.theta.00))

	# Sample from the priors for tau2inv.theta and tau2inv.gamma
	tau2inv.theta <- rgamma(n=1, shape = alpha.tau.theta, 
		rate = beta.tau.theta)
	tau2inv.gamma <- rgamma(n=1, shape = alpha.tau.gamma, 
		rate = beta.tau.gamma)
	tau.gamma <- 1/sqrt(tau2inv.gamma)
	tau.theta <- 1/sqrt(tau2inv.theta)

	# Sample sigma2inv.gamma
	sigma2inv.gamma <- rgamma(n=1, shape = alpha.sigma.gamma, 
		rate = beta.sigma.gamma)   
	sigma.gamma <- 1/sqrt(sigma2inv.gamma)
	
	# Sample the vector parameters
	mu.gamma <- rep(0, B)
	mu.theta <- rep(0, B)
	pi <- rep(0, B)
	sigma2inv.theta <- rep(0, B)
	sigma.theta <- rep(0, B)
	
	for (i in (1:B)) {
		mu.gamma[i] <- rnorm( n=1, mean = mu.gamma.0, 
			sd = 1/sqrt(tau2inv.gamma))
		mu.theta[i] <- rnorm( n=1, mean = mu.theta.0, 
			sd = 1/sqrt(tau2inv.theta))
		pi[i] <- rbeta (n=1, shape1 = alpha.pi, 
			shape2 = beta.pi)
		sigma2inv.theta[i] <- rgamma(n=1, 
			shape = alpha.sigma.theta, 
			rate = beta.sigma.theta)			
		sigma.theta[i] <- 1/sqrt(sigma2inv.theta[i])
	}

  params <- c(alpha.pi, beta.pi, mu.gamma.0, mu.theta.0,
  	tau.theta, tau.gamma, sigma.gamma,
  	mu.gamma, mu.theta, pi, sigma.theta)
	return(params)
}

generate.data.AEmodel <- function(params, inputs){
#########################################################################
##  
##  generate.data.AEmodel
##
##  given the parameters, generate data from the model
##  
##  this function is used as input to the Cook and Rubin validation
##
	
	# get the number of body systems, AEs, and treament and control 
	# patients
	B <- inputs[[1]]
	Nobs <- inputs[[2]]
	Nc <- inputs[[3]]
	Nt <- inputs[[4]]
	verbose <- inputs[[5]]
	
	sigma.gamma <- params[7]
	mu.gamma <- params[8:(B+7)]
	mu.theta <- params[(B+8):(2*B+7)]
	pi <- params[(2*B+8):(3*B+7)]
	sigma.theta <- params[(3*B+8):(4*B+7)]

	# Sample the AE-specific parameters
	maxNobs <- max(Nobs)
	c <- matrix(0, nrow = B, ncol = maxNobs)
	t <- matrix(0, nrow = B, ncol = maxNobs)
	gamma <- matrix(0, nrow = B, ncol = maxNobs)
	theta <- matrix(0, nrow = B, ncol = maxNobs)
	noDrugEffect <- matrix(0, nrow = B, ncol = maxNobs)
	theta.latent <- matrix(0, nrow = B, ncol = maxNobs)

	for (i in (1:B)) {
		for(j in (1:Nobs[i])) {
		
			noDrugEffect[i,j] <- rbinom(n=1, size=1, prob=pi[i])
			theta.latent[i,j] <- rnorm(n=1, mean = mu.theta[i], 
				sd = sigma.theta[i])
			gamma[i,j] <- rnorm(n=1, mean = mu.gamma[i], 
				sd = sigma.gamma)
			theta[i,j] <- (1-noDrugEffect[i,j]) * theta.latent[i,j] 
				
			c[i,j] <- exp(gamma[i,j]) / (1+exp(gamma[i,j]))
			logitT <- gamma[i,j] + theta[i,j]
			t[i,j] <- exp(logitT) / (1+exp(logitT))
		}
	} 

	# sample the outcome
	X <- matrix(nrow = B, ncol = maxNobs)
	Y <- matrix(nrow = B, ncol = maxNobs)
	for (i in (1:B)) {
		for(j in (1:Nobs[i])) {
			X[i, j] <- rbinom(n=1, size = Nc, prob = c[i, j])
			Y[i, j] <- rbinom(n=1, size = Nt, prob = t[i, j])
		}
	}

	aeData <- list("B"=B, 
    "Nt"=Nt, "Nc"=Nc, 
    "Nobs"=Nobs, 
    "Y"= Y,	"X" = X)

	if (F){    #(verbose){
		# print the parameter values
		print("*******************************************************************")
		print("parameter values:")
		print("")
		print("sigma.gamma:")
		print(sigma.gamma)
		print("mu.gamma")
		print(mu.gamma)
		print("mu.theta")
		print(mu.theta)
		print("pi")
		print(pi)
		print("sigma.theta")
		print(sigma.theta)
		
		# print the data
		print("*******************************************************************")
		print("data:")
		print(aeData)
	}
	
	return(aeData)
}

analyze.data.AEmodel <- function(y, params.true, inputs){
#########################################################################
##  
##  analyze.data.AEmodel
##
##  run Splus Bayes on the linear regression model to get posterior samples
##  
##  this function is used as input to the Cook and Rubin validation
##
		
	# get the hyperparameters
	hyper <- inputs[[1]]
	# get the chain specifications
	engine <-     											inputs[[2]]
	chainLength <- 											inputs[[3]]
	thinning <- 												inputs[[4]]
	burnin <- 													inputs[[5]]
	
	# get the exposures for each data point
	verbose <- 													inputs[[6]] 
	plotPoorConvergence <- 							inputs[[7]]
		
	posterior.samples <- runAEModel(aeCounts = y, nIter = chainLength,
		nThin = thinning, nBurnin = burnin, engine = engine, hyper = hyper)
	B <- y$B
	posterior.samples <- posterior.samples[,(1:(7+4*B))]

	if (verbose)	{
		# print a summary
		# print(summary(posterior.samples, maxVars = 60))
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
