{cat("----- Compare bhbm, BUGS for the seeds data w/ multiple chains ----------\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhbm
# Description: Compare the results for the bhbm function
#  to those from the BUGS engine, for a binomial-outcome
#  regression model fit to the seeds data, with multiple
#  random effects and multiple chains.


	data(seeds)
	
	# Specify the prior
	seedsPrior <- bhbm.prior ( 
	  random.var = bayes.invChisq(df=3, sigma0.sq = 1) )
	  
	inits <- list( list( random.var = rep(1, 4), fixed.coef = 0,
	    random.coef = matrix( rep(0, 8), nrow = 4 ) ),
	    list( random.var = rep(1, 4), fixed.coef = 0,
	    random.coef = matrix( rep(0, 8), nrow = 4 ) ) )
	seedsSampler <- bhbm.sampler( nSamples = 1000, nBurnin=10000,
	  nThin = 20, nChain = 2, init.point="user's choice", 
	  params.init = inits )
	  
  seeds.bhbm <- list( r = seeds$r, n = seeds$n, x1 = seeds$x1,
  	x2 = seeds$x2 )
  seeds.bhbm$group <- c( rep(1, 5), rep(2, 16) )
	
        oldopt <- options(contrasts=
		c(factor="contr.treatment", ordered="contr.poly"))
	bhbm.samples <- bhbm( data = seeds.bhbm, 
	  trials.formula = ~n, random.formula = r~x1+x2+x1:x2,
	  group.formula = ~group,
	  overdispersion = "none", prior = seedsPrior,
	  sampler = seedsSampler, debug = F )
        options(oldopt)
	
	# Prepare the data for a BUGS analysis
	seeds.bugs <- seeds.bhbm
	seeds.bugs$N <- length(seeds$n)
	
	bugsModel <- function (){
	  for( i in 1 : N ) {
	    r[i] ~ dbin(p[i],n[i])
	    logit(p[i]) <- alpha0[group[i]] + alpha1[group[i]] * x1[i] + 
	      alpha2[group[i]] * x2[i] + alpha3[group[i]] * x1[i] * x2[i]
	  }
	  for( j in 1:2 ){
		  alpha0[j] ~ dnorm(0.0,prec0)
		  alpha1[j] ~ dnorm(0.0,prec1)
		  alpha2[j] ~ dnorm(0.0,prec2)
		  alpha3[j] ~ dnorm(0.0,prec3)
		}
		prec0 ~ dgamma(1.5, 1.5)
		prec1 ~ dgamma(1.5, 1.5)
		prec2 ~ dgamma(1.5, 1.5)
		prec3 ~ dgamma(1.5, 1.5)
		sigma0 <- 1/sqrt(prec0)
		sigma1 <- 1/sqrt(prec1)
		sigma2 <- 1/sqrt(prec2)
		sigma3 <- 1/sqrt(prec3)
		invisible()
	}
		
	bugs.samples <- posteriorSamples( data = seeds.bugs,
	  parametersToSave = c( "alpha0", "alpha1", "alpha2",
	  "alpha3", "sigma0", "sigma1", "sigma2", "sigma3" ), 
	  model = bugsModel, 
	  nIter = 1000, debug = F,
	  nThin = 200, nBurnin = 30000, DIC = F, engine = "WinBUGS" )

	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index <- c(9:12, 1, 3, 5, 7, 2, 4, 6, 8)
	bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples), 
	  getSamples(bhbm.samples), print.ks=F, pvalue.cutoff = (0.01 / 100)) 

}
