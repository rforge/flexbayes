{cat("----- Compare bhpm, BUGS on the epilepsy data ----------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhpm
# Description: Compare the results for the bhpm function
#  to those from the BUGS engine, for a Poisson-outcome
#  regression model fit to the epilepsy data.

# First fit with bhpm
# Specify the prior
random.var.prior = bayes.invChisq( df = 3, sigma0.sq = 1 )
epilepsyPrior <- bhpm.prior ( 
  sigma2 = bayes.invChisq(df = 3, sigma0.sq = 1),
  random.var = random.var.prior,
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
oldopt <- options(contrasts=c(factor="contr.treatment", ordered="contr.poly"))
bhpm.samples <- bhpm( data = epilepsyBhpm, 
  random.formula = randomEff, 
  group.formula = ~ subj,
  prior = epilepsyPrior,
  sampler = epilepsySampler, overdispersion = "log-normal" )
options(oldopt)

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

  for( j in 1 : nSubj ) {
    mean[j] ~ dnorm( 0, tau2inv)
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
  parametersToSave = c( "sigma", "tau", "mean",
  "lambda" ), model = bugsModel, 
  nIter = 1000,
  nThin = 10, DIC = F, engine = "WinBUGS" )

# compare the sample distributions.
# The comparison passes for some of the lambda parameters
# but not all, because of a mixing problem with the bhpm
# log-normal model.  Here, only compare the first several lambda 
# parameters.
bugs.index <- c(296, 297, (237:295), (1:236))
bugs.samples <- bugs.samples[,bugs.index]
compareSampleDistributions(getSamples(bugs.samples, param=(1:100)), 
  getSamples(bhpm.samples, param=(1:100)), print.ks=F, 
  pvalue.cutoff = (0.01 / 100))
}
