{cat("----- Compare bhpm, BUGS gamma-conj models with multiple groups ----------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhpm
# Description: Compare the results for the bhpm function
#  to those from the BUGS engine, for a Poisson-outcome
#  regression model fit to the pumps data as 
#  described in the S+Bayes manual.

# specify the priors
xi.prior = bayes.uniformShrinkage (1)

prior.distn = bhpm.prior ( xi = xi.prior,
  level2.coef = bayes.normal(0, 1), 
  random.var = bayes.invChisq(3, 1),
  common.glm = 2 )
  
# create the chain specifications
sampler.specs <- bhpm.sampler( nSamples=1000, nThin= 100, 
  nChains = 1, nBurnin = 10000,
  init.point = "prior", update.cov = 0 ) 
  
pumps.bhpm <- pumps
pumps.bhpm$ind <- c(rep(1,5), rep(2,5))
pumps.bhpm$e <- rep(1, 10)

# run the MCMC
oldopt <- options(contrasts=c(factor="contr.treatment", ordered="contr.poly"))
bhpm.samples <- bhpm( exposure.formula= ~e, random.formula= z~1,
  level2.formula = ~1, group.formula = ~ind,
  data=pumps.bhpm, overdispersion= "gamma-conj", prior= prior.distn, 
  sampler= sampler.specs, debug = F, random.seed = 11 ) 
options(oldopt)

# Prepare the data for a BUGS analysis
Ndata <- length(pumps.bhpm$e)
pumps.bugs <- list( e = pumps.bhpm$e, z = pumps.bhpm$z,
  ind = pumps.bhpm$ind, N = Ndata)

bugsModel <- function (){

 	# define the parameters of the gamma components.  
	# The hyperparameter xi[j] should have a uniform 
	# shrinkage prior
  xi <- z0 * (shrink) / (1-shrink)
  shrink ~ dunif(0,1)
  z0 <- 1
  
	# the second-level coefficient
  alpha0 ~ dnorm( 0, 1 )
  
  # the random effects
  for( j in 1:2 ){
    beta[j] ~ dnorm( alpha0, tau2inv )
  }
  tau2inv ~ dgamma(1.5, 1.5)
  tau <- 1/sqrt(tau2inv)
  
	# for each adverse event, the poisson intensities have a common gamma prior
  for( i in 1 : N ) {
		# the number of adverse events is distributed 
		# poisson with mean expected[i] * lambda[i]
   	z[i] ~ dpois (poissonMean[i])
		poissonMean[i] <- lambda[i] * e[i]
	  lambda[i] ~ dgamma (shape[i], invScale[i])
	  shape[i] <- xi
		invScale[i] <- shape[i] / mu[i]
		log( mu[i] ) <- beta[ ind[ i ] ]
  }
  invisible()		
}

bugs.samples <- posteriorSamples( data = pumps.bugs,
  parametersToSave = c( "xi", "alpha0", "tau", "beta",
  "lambda" ), model = bugsModel, 
  nIter = 1000, nBurnin = 100000,
  nThin = 100, DIC = F, engine = "WinBUGS" )

# Compare the sample distributions.
# Adjust for the multiple test scripts being run.
bugs.index <- c(15, 1, 14, 2:13)
bugs.samples <- bugs.samples[,bugs.index]
compareSampleDistributions(getSamples(bugs.samples), 
  getSamples(bhpm.samples), print.ks=F, 
  pvalue.cutoff = (0.01 / 100))
}
