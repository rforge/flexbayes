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
	  random.var = bayes.invChisq( 3, 1 ),
	  level2.coef = bayes.normal( 0, 1 ) )
	  
	toxo.sampler <- bhbm.sampler( nSamples = 1000,
	  nThin = 100, nBurnin = 1000 )
	
        oldopt <- options(contrasts=
		c(factor="contr.treatment", ordered="contr.poly"))
	bhbm.samples <- bhbm( data = toxo.dat, 
	  trials.formula = ~ni, random.formula = yi~1,
	  level2.formula = ~1, prior = toxo.prior,
	  group.formula = ~city,
	  sampler = toxo.sampler )
        options(oldopt)
	
	# Prepare the data for a BUGS analysis
	toxo.bugs <- list(ni =toxo.dat$ni, yi = toxo.dat$yi)
	toxo.bugs$N <- length(toxo.dat$n)
	
	bugsModel <- function (){
	  for( i in 1 : N ) {
			yi[i] ~ dbin(p[i],ni[i])
			logit(p[i]) <- beta[i]
			beta[i] ~ dnorm( alpha0, tau2inv )
		}
		alpha0 ~ dnorm(0.0,1.0)
		tau2inv ~ dgamma( 1.5, 1.5 )
		tau <- 1/sqrt(tau2inv)
		invisible()
	}
	
	bugsInits <- list( alpha0 = 0, tau2inv = 1 )
	bugsInits <- list(bugsInits)
	
	bugs.samples <- posteriorSamples( data = toxo.bugs,
	  parametersToSave = c( "alpha0", "tau", "beta" ), 
	  model = bugsModel, 
	  nIter = 1000, inits = bugsInits,
	  nThin = 50, DIC = F, engine = "WinBUGS" )
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index <- c( 1, 12, (2:11) )
	bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples), 
	  getSamples(bhbm.samples), print.ks=F, pvalue.cutoff = (0.01 / 100)) 
}
