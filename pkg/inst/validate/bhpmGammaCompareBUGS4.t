{cat("----- Compare bhpm, BUGS with a uniform shrinkage prior ----------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhpm
# Description: Compare the results for the bhpm function
#  to those from the BUGS engine, for a Poisson-outcome
#  regression model fit to the pumps data as 
#  described in the S+Bayes manual.

# specify the priors
xi.prior = bayes.uniformShrinkage (1)

prior.distn = bhpm.prior ( xi = xi.prior,
  fixed.coef = bayes.normal(0, 1),
  common.glm = 1 )
  
# create the chain specifications
sampler.specs <- bhpm.sampler( nSamples=1000, nThin= 50, 
  nChains = 1, nBurnin = 1000,
  init.point = "prior", update.cov = 0 ) 
  
pumps.bhpm <- pumps
pumps.bhpm$ind <- c(rep(1,5), rep(2,5))

# run the MCMC
oldopt <- options(contrasts=c(factor="contr.treatment", ordered="contr.poly"))
bhpm.samples <- bhpm( exposure.formula= ~e, fixed.formula= z~1,
  group.formula = ~ind,
  data=pumps.bhpm, overdispersion= "gamma-conj", prior= prior.distn, 
  sampler= sampler.specs ) 
options(oldopt)

# Prepare the data for a BUGS analysis
Ndata <- length(pumps$e)
pumps.bugs <- list( e = pumps$e, z = pumps$z,
  ind = pumps.bhpm$ind, N = Ndata)

bugsModel <- function (){

 	# define the parameters of the gamma components.  
	# The hyperparameter xi[j] should have a uniform 
	# shrinkage prior
  xi[1] <- z0 * (shrink[1]) / (1-shrink[1])
  shrink[1] ~ dunif(0,1)
  xi[2] <- z0 * (shrink[2]) / (1-shrink[2])
  shrink[2] ~ dunif(0,1)
  z0 <- 1
  
	# the random coefficients
  beta0 ~ dnorm( 0, 1 )
  
	# for each adverse event, the poisson intensities have a common gamma prior
  for( i in 1 : N ) {
		# the number of adverse events is distributed 
		# poisson with mean expected[i] * lambda[i]
   	z[i] ~ dpois (poissonMean[i])
		poissonMean[i] <- lambda[i] * e[i]
	  lambda[i] ~ dgamma (shape[i], invScale[i])
	  shape[i] <- xi[ind[i]]
		invScale[i] <- shape[i] / mu[i]
		log( mu[i] ) <- beta0
  }
  invisible()		
}

bugs.samples <- posteriorSamples( data = pumps.bugs,
  parametersToSave = c( "beta0", "xi", 
  "lambda" ), model = bugsModel, 
  nIter = 1000, nBurnin = 100000,
  nThin = 1000, DIC = F, engine = "WinBUGS" )

# Compare the sample distributions.
# Adjust for the multiple test scripts being run.
bugs.index <- c( 1, (12:13), (2:11) )
bugs.samples <- bugs.samples[,bugs.index]
compareSampleDistributions(getSamples(bugs.samples), 
  getSamples(bhpm.samples), print.ks=F, 
  pvalue.cutoff = (0.01 / 100))
}
