{cat("----- Compare bhpm, BUGS; Wishart prior, gamma-conj overdisp. -------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhpm
# Description: Compare the results for the bhpm function
#  to those from the BUGS engine, using a Wishart prior for the 
#  random effect covariance matrix and gamma-conjugate overdispersion.
#

 	# specify the priors
 	xi.prior <- bayes.uniformShrinkage(1)
	random.var.prior <- bayes.invWishart( df = 3, scale = diag( c(2,2) ) )
	
	prior.distn = bhpm.prior ( xi = xi.prior, 
		random.var =  random.var.prior,
		common.glm = 2 )
	  
	# create the chain specifications
	sampler.specs <- bhpm.sampler( nSamples=1000, nThin= 50,
	  nChains = 1, nBurnin = 1000,
	  init.point = "prior", update.cov = 0 ) 
	
	# run the MCMC
        oldopt <- options(contrasts=
		c(factor="contr.treatment", ordered="contr.poly"))
	bhpm.samples <- bhpm( exposure.formula= ~e, random.formula= z~x,
	  data=pumps, prior= prior.distn, overdispersion= "gamma-conj",
	  sampler= sampler.specs ) 
        options(oldopt)
	
	# Prepare the data for a BUGS analysis
	Ndata <- length(pumps$e)
	pumps.bugs <- list( e = pumps$e, z = pumps$z, x = pumps$x,
	  N = Ndata, beta.mean = c(0,0), tau2inv.prec = diag( c(0.5,0.5) ) )
	
	bugsModel <- function (){
	
 	  # define the parameters of the gamma components.  
	  # The hyperparameter xi[j] should have a uniform 
	  # shrinkage prior
    xi <- z0 * (shrink) / (1-shrink)
    shrink ~ dunif(0,1)
    z0 <- 1
  
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
			lambda[i] ~ dgamma (shape[i], invScale[i])
	    shape[i] <- xi
		  invScale[i] <- shape[i] / mu[i]
		  log(mu[i]) <- beta[1] + beta[2] * x[i]
	  }
	  invisible()		
	}
	
	bugs.samples <- posteriorSamples( data = pumps.bugs,
	  parametersToSave = c( "tau2", "beta", "xi", "lambda" ), 
	  model = bugsModel, nIter = 1000,
	  nThin = 100, DIC = F, engine = "WinBUGS", debug=F )
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
  bugs.index <- c( 17, 13:16, 1:12 )
  bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples), 
	  getSamples(bhpm.samples), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100))

}
