{cat("-- Compare bhpm, BUGS; multiple random eff.s and gamma-conj overdisp. ----\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhpm
# Description: Compare the results for the bhpm function
#  to those from the BUGS engine, for a Poisson-outcome
#  regression model fit to the pumps data as 
#  described in the S+FlexBayes manual.  Use two fixed 
#  effects and gamma-conjugate overdispersion.
  
	# specify the priors
	xi.prior = bayes.uniformShrinkage (1)
	fixed.effects.prior = bayes.normal( mean = c(0,0), cov = diag( c(1,1) ) )
	
	prior.distn = bhpm.prior ( xi = xi.prior,
	  fixed.coef = fixed.effects.prior, common.glm = 2 )
	  
	# create the chain specifications
	sampler.specs <- bhpm.sampler( nSamples=1000, nThin= 50, 
	  nChains = 1, nBurnin = 1000,
	  init.point = "prior", update.cov = 0 ) 
	
	# run the MCMC
        oldopt <- options(contrasts=
		c(factor="contr.treatment", ordered="contr.poly"))
	bhpm.samples <- bhpm( exposure.formula= ~e, fixed.formula= z~x,
	  data=pumps, overdispersion= "gamma-conj", prior= prior.distn, 
	  sampler= sampler.specs ) 
        options(oldopt)
	
	# Prepare the data for a BUGS analysis
	Ndata <- length(pumps$e)
	pumps.bugs <- list( e = pumps$e, z = pumps$z, x = pumps$x,
	  N = Ndata )
	
	bugsModel <- function (){
	
	 	# define the parameters of the gamma components.  
		# The hyperparameter xi[j] should have a uniform 
		# shrinkage prior
	  xi <- z0 * (shrink) / (1-shrink)
	  shrink ~ dunif(0,1)
	  z0 <- 1
	 	shape <- xi
	 	
		# the random coefficients
	  beta[1] ~ dnorm(0, 1)
	  beta[2] ~ dnorm(0, 1)
	  
		# for each adverse event, the poisson intensities have a common gamma prior
	  for( i in 1 : N ) {
			# the number of adverse events is distributed 
			# poisson with mean expected[i] * lambda[i]
	   	z[i] ~ dpois (poissonMean[i])
			poissonMean[i] <- lambda[i] * e[i]
		  lambda[i] ~ dgamma (shape, invScale[i])
			invScale[i] <- xi / mu[i]
			log( mu[i] ) <- beta[1] + beta[2] * x[i]
	  }
	  invisible()		
	}
	
	bugs.samples <- posteriorSamples( data = pumps.bugs,
	  parametersToSave = c( "beta", "xi", 
	  "lambda" ), model = bugsModel, 
	  nIter = 1000,
	  nThin = 100, DIC = F, engine = "WinBUGS" )
	
	# compare the sample distributions
	# Adjust for the multiple test scripts being run.
	bugs.index <- c(1:2, 13, 3:12)
	bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples), 
	  getSamples(bhpm.samples), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100))
}
