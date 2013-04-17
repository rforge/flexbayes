runAEModel <- function(aeCounts, nIter = 5000, nBurnin = 10000,
	nThin = 10, diffusePrior = F, piSymm = T,
	engine = "OpenBUGS", hyper = NULL){

	##################################################
	## Run the model for adverse events given in 
	## Berry and Berry (2004)
	##
	
	nBodySystems <- aeCounts$B
	numAEsPerBodySystem <- aeCounts$Nobs
	maxNumAEs <- max(numAEsPerBodySystem)

	##### Specify the Berry and Berry model

	berryModel <- function (){
		for (i in 1:B) {

			for(j in 1:Nobs[i]) {
				X[i, j] ~ dbin(c[i, j], Nc)
				Y[i, j] ~ dbin(t[i, j], Nt)
			
				logit(c[i,j] )<- gamma[i,j]
				logit(t[i,j] )<- gamma[i,j] + theta[i,j]
			
				gamma[i,j]~dnorm(mu.gamma[i], sigma2inv.gamma) 
				theta[i,j] <- (1-noDrugEffect[i,j]) * theta.latent[i,j] 
				noDrugEffect[i,j] ~ dbern(pi[i])
				theta.latent[i,j] ~ dnorm(mu.theta[i], sigma2inv.theta[i])

				} 
			mu.gamma[i] ~ dnorm( mu.gamma.0, tau2inv.gamma)
			mu.theta[i] ~ dnorm( mu.theta.0, tau2inv.theta)
			pi[i] ~ dbeta (alpha.pi, beta.pi)
			
			sigma2inv.theta[i] ~ dgamma(alpha.sigma.theta, beta.sigma.theta)
			sigma.theta[i] <- 1/sqrt(sigma2inv.theta[i])
			}
		
		## The prior for sigma2inv.gamma
			
		sigma2inv.gamma ~ dgamma(alpha.sigma.gamma, beta.sigma.gamma)   
		sigma.gamma <- 1/sqrt(sigma2inv.gamma)
		
		## The priors for tau2inv.theta and tau2inv.gamma
		
		tau2inv.theta ~ dgamma(alpha.tau.theta, beta.tau.theta)
		tau2inv.gamma ~ dgamma(alpha.tau.gamma, beta.tau.gamma)
		tau.gamma <- 1/sqrt(tau2inv.gamma)
		tau.theta <- 1/sqrt(tau2inv.theta)
		
		mu.gamma.0 ~ dnorm(mu.gamma.00, tau2inv.gamma.00)
		mu.theta.0 ~ dnorm(mu.theta.00, tau2inv.theta.00)
		
		## The priors for alpha.pi and beta.pi
		
		alpha.pi <- alphaExp + 1
		alphaExp ~ dexp(lambda.alpha)
		beta.pi <- betaExp + 1
		betaExp ~ dexp(lambda.beta)
	}

	# Specify the hyperparameters
	if (!is.null(hyper)){
		aeCounts <- c(aeCounts, hyper)
	} else {
		if (diffusePrior){
	
			aeCounts$alpha.sigma.theta <- 0.5
			aeCounts$beta.sigma.theta <- 0.01
			aeCounts$alpha.sigma.gamma <- 0.01
			aeCounts$beta.sigma.gamma <- 0.01
			aeCounts$alpha.tau.theta <- 0.01
			aeCounts$beta.tau.theta <- 0.01
			aeCounts$alpha.tau.gamma <- 0.01
			aeCounts$beta.tau.gamma <- 0.01
			
			aeCounts$mu.gamma.00 <- 0
			aeCounts$tau2inv.gamma.00 <- 0.01
			aeCounts$mu.theta.00 <- 0
			aeCounts$tau2inv.theta.00 <- 0.01
		} else {
		
			aeCounts$alpha.sigma.theta <- 3
			aeCounts$beta.sigma.theta <- 1
			aeCounts$alpha.sigma.gamma <- 3
			aeCounts$beta.sigma.gamma <- 1
			aeCounts$alpha.tau.theta <- 3
			aeCounts$beta.tau.theta <- 1
			aeCounts$alpha.tau.gamma <- 3
			aeCounts$beta.tau.gamma <- 1
			
			aeCounts$mu.gamma.00 <- 0
			aeCounts$tau2inv.gamma.00 <- 0.1
			aeCounts$mu.theta.00 <- 0
			aeCounts$tau2inv.theta.00 <- 0.1
		}
		if (piSymm){
			aeCounts$lambda.alpha <- 1
			aeCounts$lambda.beta <- 1
		} else {
			aeCounts$lambda.alpha <- 0.001
			aeCounts$lambda.beta <- 1
		}
	}
	
	# Create a vector of initial values	
	gammaInit <- array (dim = c(nBodySystems, maxNumAEs))
	aInit <- array (dim = c(nBodySystems, maxNumAEs))
	deInit <- array (dim = c(nBodySystems, maxNumAEs))

	for (bodySystem in 1:nBodySystems){
		nAEsThisBodySystem <- numAEsPerBodySystem[bodySystem]
		gammaInit [bodySystem, 1:nAEsThisBodySystem] <- 0
		aInit [bodySystem, 1:nAEsThisBodySystem] <- 0
		deInit [bodySystem, 1:nAEsThisBodySystem] <- 0
	}

	initsBerry <- list(	"gamma" = gammaInit,
		"theta.latent" = aInit, "mu.gamma" = rep(0, nBodySystems), "mu.theta" = rep(0, nBodySystems), 
		"sigma2inv.gamma" = 0.1, "sigma2inv.theta" = rep(0.1, nBodySystems), "tau2inv.theta" = 0.1, "tau2inv.gamma" = 0.1, 
		"mu.gamma.0" = 0, "mu.theta.0" = 0, "pi" = rep(0.5, nBodySystems), "alphaExp" = 1,
		"noDrugEffect" = deInit,
		"betaExp" = 1)
	
	# the initial values must be a list of length = n.chains
	initsBerry <- list(initsBerry)  	

	# create a list of parameters to be saved from the mcmc run.  We will need the
	# theta parameters in order to do inference on the risks of the AEs
	parametersBerry <- c("alpha.pi", "beta.pi", "mu.gamma.0", "mu.theta.0",
		"tau.theta", "tau.gamma", "sigma.gamma", "mu.gamma", "mu.theta", "pi", "sigma.theta",
		"t", "c", "theta.latent")

	aeSim <- posteriorSamples (data = aeCounts,
		parametersToSave = parametersBerry, 
		model = berryModel, 	inits = initsBerry,
		nChains = 1, nIter = nIter, nBurnin = nBurnin, nThin = nThin, DIC = F,
		debug = F, engine = engine)
	
	aeSim
}

runEffToxCRModel <- function(	data, coefPriors, 
  debug = F){

	########################################################
	## Run the continuation ratio model for 
	## efficacy and toxicity given in Thall and Cook (2004),
	## p. 686
	## 

  # fit the toxicity and conditional efficacy models separately,
  # since they are modeled independently

  # Create a data set containing just toxicity and
  # dose data for the subjects
  nData <- length(data$outcome)
  if (length(data$patientDose) != nData)
    stop ("data inconsistent")
  toxModelData <- list(patientDoseLevel = rep(0, nData), 
    tox = as.integer(data$outcome == 2),
    trials = rep(1, nData) )
  for (i in (1:nData)){
    toxModelData$patientDoseLevel[i] <- 
      data$doseLevels[data$patientDose[i]]
  }
  
  # Create a data set containing just efficacy and 
  # dose data,
  # restricting to subjects with no toxicity 
  nonToxDose <- data$patientDose[data$outcome != 2]
  nonToxOutcome <- data$outcome[data$outcome != 2]
  nNonTox <- length( nonToxOutcome )
  effModelData <- list( 
    patientDoseLevel = rep(0, nNonTox),
    eff = as.integer(nonToxOutcome == 1),
    trials = rep(1, nNonTox) )
  for (i in (1:nNonTox)){
    effModelData$patientDoseLevel[i] <- 
      data$doseLevels[nonToxDose[i]]
  }

  # Assign independent normal priors for the fixed effect coefficients
  toxPrior <- bhbm.prior( fixed.coef = bayes.normal(
    mean = coefPriors$toxCoefMean, 
    cov = diag( coefPriors$toxCoefSD^2 ) ) )

  # Specify the priors 
  effPrior <- bhbm.prior( fixed.coef = bayes.normal(
    mean = coefPriors$effCoefMean, 
    cov = diag( coefPriors$effCoefSD^2 ) ) )
     
  # Specify the sampler control parameters
  sampler <- bhbm.sampler( nBurnin = 1000, 
    nThin = 10, nSamples = 2000, nChains = 1 )

  # Fit the toxicity and efficacy-given-no-toxicity
  # models
  postTox <- bhbm( 
    fixed.formula = tox ~ patientDoseLevel,
    trials.formula = ~ trials, data = toxModelData,
    prior = toxPrior, sampler = sampler )
  
  postEff <- bhbm( 
    fixed.formula = eff ~ patientDoseLevel,
    trials.formula = ~ trials, data = effModelData,
    prior = effPrior, sampler = sampler )

  # check the convergence of the chains
  poorlyConvergedParams <- 
    c(poorconvergenceParameters (x = postTox), 
      poorconvergenceParameters (x = postEff) )
  if (length(poorlyConvergedParams > 0)){
		warning(paste("Poor convergence of the following parameters: ", 
		  poorlyConvergedParams, sep = ""))
	}

  # Get the posterior samples from the models
  betaTox <- getSamples( postTox )
  betaEff <- getSamples( postEff )
      
  # Merge the samples from the two models.  This
  # requires that the same number of samples were
  # drawn in both cases.  Calculate the samples
  # for the marginal probability of efficacy.
  nDoses <- length(data$doseLevels)
  nDraws <- niter(postTox)
  if (nDraws != niter(postEff))
  	stop( "The number of samples must be the same in the efficacy and toxicity model runs")
  ptoxDoseSamp <- matrix(nrow = nDraws, ncol = nDoses)
  colnames(ptoxDoseSamp) <- 
    paste("ptoxDose[", (1:nDoses), "]", sep="")
  peffDoseSamp <- matrix(nrow = nDraws, ncol = nDoses)
  colnames(peffDoseSamp) <- 
    paste("peffDose[", (1:nDoses), "]", sep="")
  
  for (j in (1:nDoses)){
  
    # calculate the probability of toxicity at dose j
	  ptoxDoseSamp[,j] <- 
	    betaTox[,1] + betaTox[,2] * data$doseLevels[j]
	  ptoxDoseSamp[,j] <- 
	    exp( ptoxDoseSamp[,j] ) / ( 1 + exp( ptoxDoseSamp[,j] ) )
  
    # calculate the probability of efficacy at dose j
	  peffDoseSamp[,j] <- 
	    betaEff[,1] + betaEff[,2] * data$doseLevels[j]
	  peffDoseSamp[,j] <- 
	    exp( peffDoseSamp[,j] ) / ( 1 + exp( peffDoseSamp[,j] ) )
    peffDoseSamp[,j] <- peffDoseSamp[,j] * (1-ptoxDoseSamp[,j])
  }
  postSamples <- cbind(ptoxDoseSamp, peffDoseSamp)
    
  return(posterior(postSamples))
}

runEffToxPOModel <- function(	data, coefPriors, 
  engine = "OpenBUGS", debug = F){

	########################################################
	## Run the proportional odds model for 
	## efficacy and toxicity given in Thall and Russell (1998)
	## 

  # Create the data set for BUGS
  nData <- length( data$outcome )
  POdata <- list( patientDose = data$patientDose,
    doseLevels = data$doseLevels,
    tox = as.integer(data$outcome == 2),
    effOrTox = as.integer(data$outcome >= 1),
    nData = nData )
  POdata <- c( POdata, coefPriors )

  # Specify the model in BUGS format
  bugsModel <- function() {
    for (i in 1:nData){
      tox[i] ~ dbern( pTox[i] )
      effOrTox[i] ~ dbern( pEffOrTox[i] )
      logit(pTox[i]) <- mu + beta * doseLevels[patientDose[i]]
      logit(pEffOrTox[i]) <- mu + alpha + beta * doseLevels[patientDose[i]]
    }
    mu ~ dunif( muLower, muUpper )
    alpha ~ dunif( alphaLower, alphaUpper )
    beta ~ dunif( betaLower, betaUpper )
    invisible()
  }
  
  postEffTox <- posteriorSamples( data = POdata,
    parametersToSave = c("mu", "alpha", "beta"),
    model = bugsModel, nIter = 1000, nThin = 10,
    engine = engine )

  # check the convergence of the chains
  poorlyConvergedParams <- 
    poorconvergenceParameters (x = postEffTox)
  if (length(poorlyConvergedParams > 0)){
		warning(paste("Poor convergence of the following parameters: ", 
		  poorlyConvergedParams, sep = ""))
	}

  # Get the posterior samples from the models
  betaSamples <- getSamples( postEffTox, "beta" )
  alphaSamples <- getSamples( postEffTox, "alpha" )
  muSamples <- getSamples( postEffTox, "mu" )
      
  # Calculate the samples
  # for the marginal probability of efficacy and 
  # that of toxicity
  nDoses <- length(data$doseLevels)
  nDraws <- niter(postEffTox)
  ptoxDoseSamp <- matrix(nrow = nDraws, ncol = nDoses)
  colnames(ptoxDoseSamp) <- 
    paste("ptoxDose[", (1:nDoses), "]", sep="")
  peffDoseSamp <- matrix(nrow = nDraws, ncol = nDoses)
  colnames(peffDoseSamp) <- 
    paste("peffDose[", (1:nDoses), "]", sep="")
  
  for (j in (1:nDoses)){
  
    # calculate the probability of toxicity at dose j
	  ptoxDoseSamp[,j] <- 
	    muSamples + betaSamples * data$doseLevels[j]
	  ptoxDoseSamp[,j] <- 
	    exp( ptoxDoseSamp[,j] ) / ( 1 + exp( ptoxDoseSamp[,j] ) )
  
    # calculate the probability of efficacy at dose j
	  peffDoseSamp[,j] <- 
	    muSamples + alphaSamples + betaSamples * data$doseLevels[j]
	  peffDoseSamp[,j] <- 
	    exp( peffDoseSamp[,j] ) / ( 1 + exp( peffDoseSamp[,j] ) )
  }
  postSamples <- cbind(ptoxDoseSamp, peffDoseSamp)
    
  return(posterior(postSamples))
}

chooseNextDoseTR <- function(patientDose, outcome, 
  doseValues, doseTrans, currentDoseInd,
  lowestEff, highestTox, pEff, pTox, coefPriors,
  model = "CR", plot = T){

	####################################################
	## Choose the next dose in an adaptive drug trial;
	## here we use the method suggested by Thall and Russell 
	## (1998)
	##

	# Run the model, either the 
	# continuation ratio model or the proportional odds model
  data <- list (patientDose = patientDose, 
  	doseLevels = doseTrans, outcome = outcome)
    
  if (model == "CR")
    postSamples <- runEffToxCRModel(data = data,
      coefPriors = coefPriors)
  else if (model == "PO")
    postSamples <- runEffToxPOModel(data = data,
      coefPriors = coefPriors)
  else 
    stop("Model must be one of 'CR' (continuation ratio) or 'PO' (proportional odds)")

  # evaluate each of the possible doses
  # in terms of probability of efficacy and 
  # probability of toxicity.	
  nDose <- length(doseTrans)
  probHighEff <- rep (-1, nDose)
  probLowTox <- rep (-1, nDose)
  acceptable <- rep (F, nDose)
  piToxEst <- rep (-1, nDose)
  piEffEst <- rep (-1, nDose)
  piToxInterval <- matrix (nrow = nDose, ncol = 2)
  piEffInterval <- matrix (nrow = nDose, ncol = 2)
  
  for (i in (1:nDose)){
      
    # Obtain posterior samples for the prob of 
    # toxicity
    toxParamName <- paste("ptoxDose[", i, "]", sep = "")
    piToxSamp <- 
      getSamples(postSamples, toxParamName)
		
    # Obtain posterior samples for the prob of 
    # efficacy
    effParamName <- paste("peffDose[", i, "]", sep = "")
    piEffSamp <- 
      getSamples(postSamples, effParamName)
    
    # what is the posterior probability that the 
    # prob of efficacy is high enough?
    probHighEff[i] <- mean(piEffSamp > lowestEff)
		
    # what is the posterior probability that the prob
    # of toxicity is low enough?
    probLowTox[i] <- mean(piToxSamp < highestTox)
		
    # calculate whether or not this dose is 
    # acceptable
    acceptable[i] <- (probHighEff[i] > pEff) & 
      (probLowTox[i] > pTox)
      
    # calculate point and interval estimates of the 
    # prob of efficacy and toxicity
    piToxEst[i] <- mean(piToxSamp)
    piEffEst[i] <- mean(piEffSamp)
    piToxInterval[i,] <- quantile(x = piToxSamp, probs = c(0.1, 0.9))
    piEffInterval[i,] <- quantile(x = piEffSamp, probs = c(0.1, 0.9))
  }
	
  if (sum(acceptable) == 0){
    chosenDoseInd <- -1
  } else {
    # find the estimated best dose out of the 
    # acceptable doses.  If there is 
    # more than one best dose, choose the lowest
    # one.
    bestDoseEst <- which(
      probHighEff == max(probHighEff * acceptable))
    bestDoseEst <- min(bestDoseEst)

    # if the estimated best dose is more than one 
    # dose level above the current level, increment
    # the current level by one.  Otherwise, choose
    # the estimated best dose.
    if (bestDoseEst > (currentDoseInd + 1)){
      chosenDoseInd <- currentDoseInd + 1
    } else {
      chosenDoseInd <- bestDoseEst
    }
  }
  # create a plot of the intervals for piEff and piTox
  if (plot) {
    effToxIntervalPlot (piEffEst = piEffEst, 
      piToxEst = piToxEst, piEffInterval = piEffInterval,
      piToxInterval = piToxInterval, 
      lowestEff = lowestEff, highestTox = highestTox,
      doseValues = doseValues)
  }
    
  return(chosenDoseInd)
}

effToxIntervalPlot <- function(piEffEst, piToxEst,
	piEffInterval, piToxInterval,
	lowestEff, highestTox, doseValues){
	
	################################################
	##
	## Create a plot of the point and interval 
	## estimates for toxicity and efficacy 
	## at the different dose levels
	##
	
	# Turn off warnings during the plotting; this 
	# avoids the printing of an unnecessary warning.
	oldWarn <- options("warn")$warn
	options(warn=-1)
	on.exit(options(warn=oldWarn))
	
	par(mfrow=c(2,1))
	error.bar(x=doseValues, y=piEffEst, 
		lower = piEffInterval[,1], 
		upper = piEffInterval[,2], incr=F, 
		ylim=c(0,1), xaxt='n',
		ylab='Pr(Efficacy)',xlab='Dose',
		main='Pr(Efficacy) by Dose')
		
	axis(1, at=doseValues, labels=doseValues)
	abline(h=lowestEff,lty=2) 

	key( text = "Point and 80-percent interval estimates",
	  corner = c(0,1), x = min( doseValues ),
	  y = 1.0 )

	error.bar(x=doseValues, y=piToxEst, 
		lower = piToxInterval[,1], 
		upper = piToxInterval[,2], incr=F, 
		ylim=c(0,1), xaxt='n',
		ylab='Pr(Toxicity)',xlab='Dose',
		main='Pr(Toxicity) by Dose')
	axis(1, at=doseValues, labels=doseValues)
	abline(h=highestTox,lty=2)
	
	key( text = "Point and 80-percent interval estimates",
	  corner = c(0,1), x = min( doseValues ),
	  y = 1.0 )
	
}
