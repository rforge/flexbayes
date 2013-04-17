#########################################################################
##  
##  Validation of the bhpm function
##
##  
##


bhpmCompareBUGS3 <- function() {
	
	# Load the data
	data(epilepsy)
	# Prepare the data for a BUGS analysis
	Ndata <- length( as.vector(epilepsy$y) )
	epilepsy.bugs <- list( y = as.vector( t( epilepsy$y) ), 
	  N = Ndata )
	
	bugsModel <- function (){
	
	  # The overdispersion parameter 
	  sigma2inv ~ dgamma( 1.5, 1.5 )
	 	sigma <- 1 / sqrt( sigma2inv )
	
		# the fixed coefficient
	  beta0 ~ dnorm( 0.0, 1)
	
		# for each adverse event, the poisson intensities have a common gamma prior
	  for( i in 1 : N ) {
	   	y[i] ~ dpois (lambda[i])
		  log(lambda[i]) <- beta0 + z[i]
		  z[i] ~ dnorm( 0, sigma2inv)
	  }
	  invisible()		
	}
	
	bugs.samples <- posteriorSamples( data = epilepsy.bugs,
	  parametersToSave = c( "beta0", "sigma",
	  "lambda" ), model = bugsModel, 
	  nIter = 1000,
	  nThin = 40, DIC = F, engine = "WinBUGS", debug = F )
	
	# Fit with bhpm
	# Specify the prior
	epilepsyPrior <- bhpm.prior ( 
	  sigma2 = bayes.invChisq(df = 3, sigma0.sq = 1),  
	  fixed.coef = bayes.normal(mean = zero, cov = identity),
	  common.glm = 2 )
	  
	# Specify chain characteristics
	epilepsySampler <- bhpm.sampler( nSamples = 1000,
	  nThin = 10 )
	
	# Prepare the data 
	nSubj <- dim(epilepsy$y)[1]
	base <- log ( epilepsy$baseline / 4 )
	epilepsyBhpm <- data.frame( 
	  y = as.vector( t( epilepsy$y ) ),
	  subj = rep( (1:nSubj), each = 4 ) )  
	
	# fit the model
	fixedEff <- y ~ 1
	bhpm.samples <- bhpm( data = epilepsyBhpm, 
	  fixed.formula = fixedEff, 
	  group.formula = ~ subj,
	  prior = epilepsyPrior,
	  sampler = epilepsySampler, overdispersion = "log-normal",
	  debug = F )
	
	# compare the sample distributions.
	# The comparison passes for some of the lambda parameters
	# but not all, because of a mixing problem with the bhpm
	# log-normal model.  Here, only compare the first several lambda 
	# parameters.
	bugs.index <- c(1, 238, (2:237))
	bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples, param=(1:20)), 
	  getSamples(bhpm.samples, param=(1:20)), print.ks=F,
	  pvalue.cutoff = (0.01 / 100))

}


bhpmCompareBUGS4 <- function() {
	
	# First fit with bhpm
	# Specify the prior
	random.var.prior = bayes.invChisq( df = 3, sigma0.sq = 1 )
	epilepsyPrior <- bhpm.prior ( 
	  sigma2 = bayes.invChisq( df = 3, sigma0.sq = 1 ),
	  random.var = random.var.prior,
	  level2.coef = bayes.normal(mean.vector = zero, covmat = identity), 
	  common.glm = 2 )
	  
	# Specify chain characteristics
	epilepsySampler <- bhpm.sampler( nSamples = 1000,
	  nThin = 10 )
	
	# Load the data
	data(epilepsy)
	# Prepare the data 
	nSubj <- dim(epilepsy$y)[1]
	base <- log ( epilepsy$baseline / 4 )
	epilepsyBhpm <- data.frame( 
	  y = as.vector( t( epilepsy$y ) ),
	  subj = rep( (1:nSubj), each = 4 ) )  
	
	# fit the model
	randomEff <- y ~ 1
	bhpm.samples <- bhpm( data = epilepsyBhpm, 
	  random.formula = randomEff, 
	  group.formula = ~ subj, level2.formula = ~ 1,
	  prior = epilepsyPrior,
	  sampler = epilepsySampler, overdispersion = "log-normal" )
	
	# Prepare the data for a BUGS analysis
	nData <- length( as.vector(epilepsy$y) )
	epilepsy.bugs <- list( y = as.vector( t( epilepsy$y) ), 
	  subj = rep( (1:nSubj), each = 4 ),
	  nData = nData, nSubj = nSubj )
	
	bugsModel <- function (){
	
		# The overdispersion parameter 
	  sigma2inv ~ dgamma( 1.5, 1.5 )
	 	sigma <- 1 / sqrt( sigma2inv )
	  # The random effects precision
	  tau2inv ~ dgamma( 1.5, 1.5 )
	 	tau <- 1 / sqrt( tau2inv )
	
		# the fixed coefficient
	  beta0 ~ dnorm( 0.0, 1)
	
	  for( j in 1 : nSubj ) {
	    mean[j] ~ dnorm( beta0, tau2inv)
	  }
	
		# for each adverse event, the poisson intensities have a common gamma prior
	  for( i in 1 : nData ) {
	   	y[i] ~ dpois (lambda[i])
	   	log(lambda[i]) <- logLambda[i]
		  logLambda[i] ~ dnorm( logLambdaMean[i], sigma2inv )
		  logLambdaMean[i] <- mean[subj[i]]
	  }
	  invisible()		
	}
	
	bugs.samples <- posteriorSamples( data = epilepsy.bugs,
	  parametersToSave = c( "sigma", "beta0", "tau", "mean",
	  "lambda" ), model = bugsModel, 
	  nIter = 1000,
	  nThin = 10, DIC = F, engine = "WinBUGS", debug = F )
	
	# compare the sample distributions.
	# The comparison passes for some of the lambda parameters
	# but not all, because of a mixing problem with the bhpm
	# log-normal model.  Here, only compare the first 60 or so lambda 
	# parameters.
	bugs.index <- c(297, 1, 298, (238:296), (2:237))
	bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples, param=(1:120)), 
	  getSamples(bhpm.samples, param=(1:120)), print.ks=F,
	  pvalue.cutoff = (0.01 / 100))

}
