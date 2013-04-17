{cat("----- Compare bhpm, BUGS on the pumps data ----------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhpm
# Description: Compare the results for the bhpm function
#  to those from the BUGS engine, for a Poisson-outcome
#  regression model fit to the pumps data as 
#  described in the S+FlexBayes manual.

# specify the priors
prior.distn = bhpm.prior( fixed.coef = bayes.normal( 0, 1 ) )
  
# create the chain specifications
sampler.specs <- bhpm.sampler( nSamples=1000, nThin= 100, 
  nChains = 1, nBurnin = 1000,
  init.point = "prior", update.cov = 0 ) 

# run the MCMC
oldopt <- options(contrasts=c(factor="contr.treatment", ordered="contr.poly"))
bhpm.samples <- bhpm( exposure.formula= ~e, fixed.formula= z~1,
  data=pumps, prior= prior.distn, 
  sampler= sampler.specs, debug = F ) 
options(oldopt)

# Prepare the data for a BUGS analysis
Ndata <- length(pumps$e)
pumps.bugs <- list( e = pumps$e, z = pumps$z,
  N = Ndata )

bugsModel <- function (){
	# the random coefficients
  beta0 ~ dnorm( 0.0, 1.0 )

	# for each adverse event, the poisson intensities have a common gamma prior
  for( i in 1 : N ) {
   	z[i] ~ dpois (poissonMean[i])
		poissonMean[i] <- lambda[i] * e[i]
	  log( lambda[i] ) <- beta0
  }
  invisible()		
}

bugs.samples <- posteriorSamples( data = pumps.bugs,
  parametersToSave = c( "beta0" ), model = bugsModel, 
  nIter = 1000,
  nThin = 100, DIC = F, engine = "WinBUGS" )

# Compare the sample distributions.
# Adjust for the multiple test scripts being run.
compareSampleDistributions(getSamples(bugs.samples), 
  getSamples(bhpm.samples), print.ks=F, 
  pvalue.cutoff = (0.01 / 100))
}
