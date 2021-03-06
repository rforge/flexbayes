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
sampler.specs <- bhpm.sampler( nSamples = 1000, nThin= 100, 
  nChains = 1, nBurnin = 100,
  init.point = "prior", update.cov = 0 ) 

pumps.bhpm <- pumps
pumps.bhpm$ind <- (1:length(pumps$x))

# run the MCMC
oldopt <- options(contrasts=c(factor="contr.treatment", ordered="contr.poly"))
bhpm.samples <- bhpm( exposure.formula= ~e, random.formula= z~x,
  level2.formula = ~1, group.formula = ~ind, data=pumps.bhpm, 
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
  tau2inv0 ~ dgamma( 1.5, 1.5 )
  tau0 <- 1/sqrt(tau2inv0)
  tau2inv1 ~ dgamma( 1.5, 1.5 )
  tau1 <- 1/sqrt(tau2inv1)

	# for each adverse event, the poisson intensities have a common gamma prior
  for( i in 1 : N ) {
   	z[i] ~ dpois (poissonMean[i])
		poissonMean[i] <- lambda[i] * e[i]
	  log( lambda[i] ) <- beta0[i] + beta1[i] * x[i]
	  beta0[i] ~ dnorm( alpha0, tau2inv0 )
	  beta1[i] ~ dnorm( alpha1, tau2inv1 )
  }
  invisible()		
}

bugs.samples <- posteriorSamples( data = pumps.bugs,
  parametersToSave = c( "alpha0", "alpha1", "beta0", 
  "beta1", "tau0", "tau1" ), model = bugsModel, nIter = 1000,
  nThin = 500, DIC = F, engine = "WinBUGS" )

# Compare the sample distributions.
# Adjust for the multiple test scripts being run.
bugs.index <- c( (1:2), 23:24, 3, 13, 4, 14, 5, 15, 6, 16,
  7, 17, 8, 18, 9, 19, 10, 20, 11, 21, 12, 22 )
bugs.samples <- bugs.samples[,bugs.index]
compareSampleDistributions(getSamples(bugs.samples), 
  getSamples(bhpm.samples), print.ks=F, 
  pvalue.cutoff = (0.01 / 100))
}
