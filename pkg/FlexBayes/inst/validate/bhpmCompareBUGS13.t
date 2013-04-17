{cat("----- Compare bhpm, BUGS; Wishart prior for rand. covar. -------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhpm
# Description: Compare the results for the bhpm function
#  to those from the BUGS engine, using a Wishart prior for the 
#  random effect covariance matrix.
#

 	# specify the priors
	random.var.prior = bayes.invWishart( df = 3, scale = diag( c(2,2) ) )
	
	prior.distn = bhpm.prior ( random.var =  random.var.prior )
	  
	# create the chain specifications
	sampler.specs <- bhpm.sampler( nSamples=100, nThin= 50,
	  nChains = 1, nBurnin = 1,
	  init.point = "prior", update.cov = 0 ) 
	
	# run the MCMC
        oldopt <- options(contrasts=
		c(factor="contr.treatment", ordered="contr.poly"))
	bhpm.samples <- bhpm( exposure.formula= ~e, random.formula= z~x,
	  data=pumps, prior= prior.distn, 
	  sampler= sampler.specs ) 
        options(oldopt)
	
	# Prepare the data for a BUGS analysis
	Ndata <- length(pumps$e)
	pumps.bugs <- list( e = pumps$e, z = pumps$z, x = pumps$x,
	  N = Ndata, beta.mean = c(0,0), tau2inv.prec = diag( c(0.5,0.5) ) )
	
	bugsModel <- function (){
	
	 	tau2inv[1:2,1:2] ~ dwish( tau2inv.prec[,], 3 )
	 	tau2[1:2,1:2] <- inverse( tau2inv[,] )
	
		# the random coefficients
	  beta[1:2] ~ dmnorm( beta.mean[], tau2inv[,])
	  
		# for each adverse event, the poisson intensities have a common gamma prior
	  for( i in 1 : N ) {
			# the number of adverse events is distributed 
			# poisson with mean expected[i] * lambda[i]
	   	z[i] ~ dpois (poissonMean[i])
			poissonMean[i] <- lambda[i] * e[i]
		  log( lambda[i] ) <- beta[1] + beta[2] * x[i]
	  }
	  invisible()		
	}
	
	bugs.samples <- posteriorSamples( data = pumps.bugs,
	  parametersToSave = c( "tau2", "beta" ), model = bugsModel, 
	  nIter = 1000,
	  nThin = 100, DIC = F, engine = "WinBUGS" )
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
  bugs.index <- c( 3:6, 1:2 )
  bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples), 
	  getSamples(bhpm.samples), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100))
}
