{cat("----- Compare bhbm, BUGS for the toxo data w/o overdispersion ----------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhbm
# Description: Compare the results for the bhbm function
#  to those from the BUGS engine, for a binomial-outcome
#  regression model fit to the toxo data, with no
#  overdispersion.
#

	# Specify the prior
	toxo.prior <- bhbm.prior( 
	  fixed.coef = bayes.normal( 0, 1 ) )
	  
	toxo.sampler <- bhbm.sampler( nSamples = 1000,
	  nThin = 100, nBurnin = 1000 )
	
        oldopt <- options(contrasts=
		c(factor="contr.treatment", ordered="contr.poly"))
	bhbm.samples <- bhbm( data = toxo.dat, 
	  trials.formula = ~ni, fixed.formula = yi~1,
	  prior = toxo.prior,
	  sampler = toxo.sampler )
        options(oldopt)
	
	# Prepare the data for a BUGS analysis
	toxo.bugs <- list(ni =toxo.dat$ni, yi = toxo.dat$yi)
	toxo.bugs$N <- length(toxo.dat$n)
	
	bugsModel <- function (){
	  for( i in 1 : N ) {
			yi[i] ~ dbin(p[i],ni[i])
			logit(p[i]) <- alpha0
		}
		alpha0 ~ dnorm(0.0,1.0)
		invisible()
	}
	
	bugsInits <- list(alpha0 = 0)
	bugsInits <- list(bugsInits)
	
	bugs.samples <- posteriorSamples( data = toxo.bugs,
	  parametersToSave = c( "alpha0" ), 
	  model = bugsModel, 
	  nIter = 1000, inits = bugsInits,
	  nThin = 50, DIC = F, engine = "WinBUGS" )
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	compareSampleDistributions(getSamples(bugs.samples), 
	  getSamples(bhbm.samples), print.ks=F, pvalue.cutoff = (0.01 / 100)) 
}
