{cat("----- Compare bhpm, BUGS on the pumps data ----------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhpm
# Description: Compare the results for the bhpm function
#  to those from the BUGS engine, for a Poisson-outcome
#  regression model fit to the pumps data as 
#  described in the S+FlexBayes manual.

# specify the priors
prior.distn = bhpm.prior( random.var = bayes.invChisq( 3, 1 ),
  level2.coef = bayes.normal( rep(0, 2), diag(2) ) )
  
# create the chain specifications
sampler.specs <- bhpm.sampler( nSamples=1000, nThin= 100, 
  nChains = 1, nBurnin = 100,
  init.point = "prior", update.cov = 0 ) 

pumps.bhpm <- pumps
pumps.bhpm$ind <- (1:length(pumps$x))

# run the MCMC
oldopt <- options(contrasts=c(factor="contr.treatment", ordered="contr.poly"))
bhpm.samples <- bhpm( exposure.formula= ~e, random.formula= z~1,
  level2.formula = ~x, group.formula = ~ind, data=pumps.bhpm, 
  prior= prior.distn, 
  sampler= sampler.specs, debug = F ) 
options(oldopt)

# Prepare the data for a BUGS analysis
Ndata <- length(pumps$e)
pumps.bugs <- list( e = pumps$e, z = pumps$z,
  x = pumps$x, N = Ndata )

bugsModel <- function (){
	# the regression coefficients
  alpha0 ~ dnorm( 0.0, 1.0 )
  alpha1 ~ dnorm( 0.0, 1.0 )
  
  # the random effects
  tau2inv ~ dgamma( 1.5, 1.5 )
  tau <- 1/sqrt(tau2inv)

	# for each adverse event, the poisson intensities have a common gamma prior
  for( i in 1 : N ) {
   	z[i] ~ dpois (poissonMean[i])
		poissonMean[i] <- lambda[i] * e[i]
	  log( lambda[i] ) <- beta[i]
	  beta[i] ~ dnorm( betaMean[i], tau2inv )
	  betaMean[i] <- alpha0 + alpha1 * x[i]
  }
  invisible()		
}

bugs.samples <- posteriorSamples( data = pumps.bugs,
  parametersToSave = c( "alpha0", "alpha1", "beta", 
  "tau" ), model = bugsModel, nIter = 1000,
  nThin = 100, DIC = F, engine = "WinBUGS" )

# Compare the sample distributions.
# Adjust for the multiple test scripts being run.
bugs.index <- c( (1:2), 13, (3:12) )
bugs.samples <- bugs.samples[,bugs.index]
compareSampleDistributions(getSamples(bugs.samples), 
  getSamples(bhpm.samples), print.ks=F, 
  pvalue.cutoff = (0.01 / 100))
}
