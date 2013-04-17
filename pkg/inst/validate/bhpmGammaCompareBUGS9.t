{cat("--- Compare bhpm, BUGS gamma-conj model w/ random eff.s, specify inits -----\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhpm
# Description: Compare the results for the bhpm function
#  to those from the BUGS engine, for a Poisson-outcome
#  regression model fit to the pumps data as 
#  described in the S+Bayes manual.  Use multiple random
#  effects and gamma-conjugate overdispersion, and 
#  specify the initial values.
#

# specify the priors
xi.prior = bayes.uniformShrinkage (1)

prior.distn = bhpm.prior ( xi = xi.prior,
  level2.coef = bayes.normal(rep(0,2), rep(1,2)), 
  random.var = bayes.invChisq(3, 1),
  common.glm = 2 )
  
# create the chain specifications
sampler.specs <- bhpm.sampler( nSamples=1000, nThin= 100, 
  nChains = 1, nBurnin = 10000,
  init.point = "user's choice", update.cov = 0,
  params.init = list( xi = 1, 
  random.coef = matrix( rep(0,4), nrow=2 ),
  level2.coef = rep(0,2), random.var = rep(1,2),
  fixed.coef = 0 )  ) 
  
pumps.bhpm <- pumps
pumps.bhpm$ind <- c(rep(1,5), rep(2,5))
pumps.bhpm$e <- rep(1, 10)

# run the MCMC
oldopt <- options(contrasts=c(factor="contr.treatment", ordered="contr.poly"))
bhpm.samples <- bhpm( exposure.formula= ~e, random.formula= z~x,
  level2.formula = ~1, group.formula = ~ind,
  data=pumps.bhpm, overdispersion= "gamma-conj", prior= prior.distn, 
  sampler= sampler.specs, debug = F, random.seed = 11 ) 
options(oldopt)

# Prepare the data for a BUGS analysis
Ndata <- length(pumps.bhpm$e)
pumps.bugs <- list( e = pumps.bhpm$e, z = pumps.bhpm$z,
  x = pumps.bhpm$x,
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
  alpha1 ~ dnorm( 0, 1 )
  
  # the random effects
  for( j in 1:2 ){
    beta0[j] ~ dnorm( alpha0, tau2inv0 )
    beta1[j] ~ dnorm( alpha1, tau2inv1 )
  }
  tau2inv0 ~ dgamma(1.5, 1.5)
  tau0 <- 1/sqrt(tau2inv0)
  tau2inv1 ~ dgamma(1.5, 1.5)
  tau1 <- 1/sqrt(tau2inv1)
  
	# for each adverse event, the poisson intensities have a common gamma prior
  for( i in 1 : N ) {
		# the number of adverse events is distributed 
		# poisson with mean expected[i] * lambda[i]
   	z[i] ~ dpois (poissonMean[i])
		poissonMean[i] <- lambda[i] * e[i]
	  lambda[i] ~ dgamma (shape[i], invScale[i])
	  shape[i] <- xi
		invScale[i] <- shape[i] / mu[i]
		log( mu[i] ) <- beta0[ ind[ i ] ] + beta1[ ind[ i ] ] * x[i]
  }
  invisible()		
}

bugs.samples <- posteriorSamples( data = pumps.bugs,
  parametersToSave = c( "xi", "alpha0", "alpha1", 
    "tau0", "tau1", "beta0", "beta1", "lambda" ), 
  model = bugsModel, 
  nIter = 1000, nBurnin = 100000,
  nThin = 100, DIC = F, engine = "WinBUGS" )

# Compare the sample distributions.
# Adjust for the multiple test scripts being run.
bugs.index <- c(19, 1:2, 17:18, 3, 5, 4, 6:16)
bugs.samples <- bugs.samples[,bugs.index]
compareSampleDistributions(getSamples(bugs.samples), 
  getSamples(bhpm.samples), print.ks=F, 
  pvalue.cutoff = (0.01 / 100))
}
