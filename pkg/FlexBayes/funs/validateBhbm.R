#########################################################################
##  
##  Validation of the bhbm function
##
##  
##


bhbmBUGScompare <- function(){
	
	data(seeds)
	
	# Specify the prior
	seedsPrior <- bhbm.prior ( 
	  sigma2 = bayes.invChisq(df = 3, sigma0.sq = 2),
	  fixed.coef = bayes.normal(mean = zero, 
	  cov = identity), common.glm = 2 )
	  
	seedsSampler <- bhbm.sampler( nSamples = 1000,
	  nThin = 10 )
	
	bhbm.samples <- bhbm( data = seeds, 
	  trials.formula = ~n, fixed.formula = r~x1+x2+x1*x2,
	  overdispersion = "logit-normal", prior = seedsPrior,
	  sampler = seedsSampler )
	
	
	# Prepare the data for a BUGS analysis
	seeds.bugs <- seeds
	seeds.bugs$N <- length(seeds$n)
	
	bugsModel <- function (){
	  for( i in 1 : N ) {
	    r[i] ~ dbin(p[i],n[i])
	    b[i] ~ dnorm(0.0,tau)
	    logit(p[i]) <- alpha0 + alpha1 * x1[i] + alpha2 * x2[i] + alpha12 * x1[i] * x2[i] + b[i]
	  }
		alpha0 ~ dnorm(0.0,1)
		alpha1 ~ dnorm(0.0,1)
		alpha2 ~ dnorm(0.0,1)
		alpha12 ~ dnorm(0.0,1)
		tau ~ dgamma(1.5, 3)
		sigma <- 1/sqrt(tau)
		invisible()
	}
	
	bugsInits <- list(alpha0 = 0, alpha1 = 0, alpha2 = 0, alpha12 = 0, tau = 1)
	bugsInits <- list(bugsInits)
	
	bugs.samples <- posteriorSamples( data = seeds.bugs,
	  parametersToSave = c( "alpha0", "alpha1", "alpha2", 
	  "alpha12", "sigma", "p" ), model = bugsModel, 
	  nIter = 1000, inits = bugsInits,
	  nThin = 50, DIC = F, engine = "WinBUGS" )

	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index <- c(1:2, 4, 3, 26, 5:25)
	bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples), 
	  getSamples(bhbm.samples), print.ks=F, pvalue.cutoff = (0.01 / 100)) 

}


bhbmBUGScompare2 <- function(){
	
	# Specify the prior
	toxo.prior <- bhbm.prior( xi = 
	  bayes.uniformShrinkage( 1 ), 
	  random.var = bayes.invChisq( df=3, sigma0.sq=1 ),
	  common.glm = 1 )
	  
	toxo.sampler <- bhbm.sampler( nSamples = 1000,
	  nThin = 100, nBurnin = 10000, init.point = "user's choice",
	  params.init = list( xi = 20, level2.coef = 0.2, 
	  fixed.coef = 0,   # this will be ignored b/c there are no fixed effects in the model
	  random.coef = 0.2, random.var = 1 ) )
	
	bhbm.samples <- bhbm( data = toxo.dat, 
	  trials.formula = ~ni, random.formula = yi~1,
	  prior = toxo.prior,
	  sampler = toxo.sampler,
	  overdispersion = "beta-conj" )
	
	
	# Prepare the data for a BUGS analysis
	toxo.bugs <- list(ni =toxo.dat$ni, yi = toxo.dat$yi)
	toxo.bugs$N <- length(toxo.dat$n)
	
	bugsModel <- function (){
	  for( i in 1 : N ) {
			yi[i] ~ dbin(p[i],ni[i])
			p[i] ~ dbeta(a[i], b[i])
			a[i] <- xi*mu[i]
			b[i] <- xi*(1-mu[i])
			logit(mu[i]) <- alpha0
		}
		alpha0 ~ dnorm(0.0,tau)
		tau ~ dgamma(1.5, 1.5)
		xi <- z0 * (shrink) / (1-shrink)
	  shrink ~ dunif(0,1)
	  z0 <- 1
		sigma <- 1/sqrt(tau)
	  invisible()
	}
	
	bugsInits <- list(alpha0 = 0, shrink = 0.5, tau = 1)
	bugsInits <- list(bugsInits)
	
	bugs.samples <- posteriorSamples( data = toxo.bugs,
	  parametersToSave = c( "xi", "sigma", "alpha0", "p" ), 
	  model = bugsModel, 
	  nIter = 1000, inits = bugsInits,
	  nThin = 50, DIC = F, engine = "WinBUGS" )
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index <- c(13, 12, 1:11)
	bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples), 
	  getSamples(bhbm.samples), print.ks=F, pvalue.cutoff = (0.01 / 100)) 

}


bhbmBUGScompare3 <- function(){
	
	data(seeds)
	
	# Specify the prior
	seedsPrior <- bhbm.prior ( 
	  sigma2 = bayes.invChisq(df = 3, sigma0.sq = 1),
	  random.var = bayes.invChisq(df = 3, sigma0.sq = 1),
	  level2.coef = bayes.normal(mean = zero, cov = identity),
	  common.glm = 2 )
	  
	seedsSampler <- bhbm.sampler( nSamples = 1000,
	  nThin = 10 )
	
	seeds.bhbm <- seeds
	seeds.bhbm$ind <- c( rep(1, 5), rep(2, 6), rep(3, 5),
	  rep(4, 5) )
	
	bhbm.samples <- bhbm( data = seeds.bhbm, 
	  trials.formula = ~n, random.formula = r~1,
	  group.formula = ~ind, level2.formula = ~x1+x2,
	  overdispersion = "logit-normal", prior = seedsPrior,
	  sampler = seedsSampler )
	
	
	# Prepare the data for a BUGS analysis
	seeds.bugs <- list( r= seeds.bhbm$r, 
	  n = seeds.bhbm$n, ind = seeds.bhbm$ind,
	  x1 = c(0,0,1,1), x2 = c(0,1,0,1) )
	seeds.bugs$N <- length(seeds$n)
	
	bugsModel <- function (){
	  for( i in 1 : N ) {
	    r[i] ~ dbin(p[i],n[i])
	    b[i] ~ dnorm(beta[ ind[ i ] ] ,sigma2inv)
	    logit(p[i]) <- b[i]
	  }
	  for( j in 1 : 4 ){
	    beta[j] ~ dnorm( betaMean[j], tau2inv )
	    betaMean[j] <-   alpha0 + alpha1 * x1[j] + alpha2 * x2[j]
	  }
	
		alpha0 ~ dnorm(0.0,1)
		alpha1 ~ dnorm(0.0,1)
		alpha2 ~ dnorm(0.0,1)
		sigma2inv ~ dgamma(1.5, 1.5)
		sigma <- 1/sqrt(sigma2inv)
		tau2inv ~ dgamma(1.5, 1.5)
		tau <- 1/sqrt(tau2inv)
		invisible()
	}
	
	bugsInits <- list(alpha0 = 0, alpha1 = 0, alpha2 = 0, 
	  tau2inv = 1, sigma2inv = 1 )
	bugsInits <- list(bugsInits)
	
	bugs.samples <- posteriorSamples( data = seeds.bugs,
	  parametersToSave = c( "sigma", "alpha0", "alpha1", "alpha2", 
	  "tau", "beta", "p" ), model = bugsModel, 
	  nIter = 1000, inits = bugsInits,
	  nThin = 50, DIC = F, engine = "WinBUGS" )
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index <- c( 29, (1:3), 30, (4:7), (8:28) )
	bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples), 
	  getSamples(bhbm.samples), print.ks=F, pvalue.cutoff = (0.01 / 100)) 

}

bhbmBUGScompare4 <- function(){
	
	data(seeds)
	
	# Specify the prior
	seedsPrior <- bhbm.prior ( 
	  random.var = bayes.invChisq(df=3, sigma0.sq = 1) )
	  
	seedsSampler <- bhbm.sampler( nSamples = 1000, nBurnin=10000,
	  nThin = 20 )
	  
  seeds.bhbm <- list( r = seeds$r, n = seeds$n, x1 = seeds$x1 )
  seeds.bhbm$group <- c( rep(1, 5), rep(2, 16) )
	
	bhbm.samples <- bhbm( data = seeds.bhbm, 
	  trials.formula = ~n, random.formula = r~x1,
	  group.formula = ~group,
	  overdispersion = "none", prior = seedsPrior,
	  sampler = seedsSampler, debug = F )
	
	
	# Prepare the data for a BUGS analysis
	seeds.bugs <- seeds.bhbm
	seeds.bugs$N <- length(seeds$n)
	
	bugsModel <- function (){
	  for( i in 1 : N ) {
	    r[i] ~ dbin(p[i],n[i])
	    logit(p[i]) <- alpha0[group[i]] + alpha1[group[i]] * x1[i]
	  }
	  for( j in 1:2 ){
		  alpha0[j] ~ dnorm(0.0,prec0)
		  alpha1[j] ~ dnorm(0.0,prec1)
		}
		prec0 ~ dgamma(1.5, 1.5)
		prec1 ~ dgamma(1.5, 1.5)
		sigma0 <- 1/sqrt(prec0)
		sigma1 <- 1/sqrt(prec1)
		invisible()
	}
		
	bugs.samples <- posteriorSamples( data = seeds.bugs,
	  parametersToSave = c( "alpha0", "alpha1",
	  "sigma0", "sigma1" ), model = bugsModel, 
	  nIter = 1000, debug = F,
	  nThin = 200, nBurnin = 30000, DIC = F, engine = "WinBUGS" )

	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index <- c(5:6, 1, 3, 2, 4)
	bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples), 
	  getSamples(bhbm.samples), print.ks=F, pvalue.cutoff = (0.01 / 100)) 

}


bhbmBUGScompareSalamanders <- function(){
	
	# Specify the prior
	salPrior <- bhbm.prior ( 
	  sigma2 = bayes.invChisq(df = 3, sigma0.sq = 1),
	  fixed.coef = bayes.normal( mean = rep(0, 4), 
	  cov = rep(1, 4) ),
	  random.var = bayes.invChisq(df = 3, sigma0.sq=1),
	  common.glm = 2 )
	
	data(salamanders)
	
	salSampler <- bhbm.sampler( nSamples = 1000,
	  nThin = 1) #50 )
	
	# Load the data
	data(salamanders)
	# Prepare the data 
	# Just analyze the first experiment
	nMates <- dim(salamanders$y)[1]
	salBhbm <- data.frame( y = salamanders$y[,1],
	  x = salamanders$x, z = salamanders$z,
	  trials = rep( 1, nMates ) )
	
	# fit the model
	fixedEff <- 
	  y ~ x.R.R + x.R.W + x.W.R + x.W.W - 1
	randomEff <- ~ z.fR1 + z.fR2 + z.fR3 + z.fR4 + 
	  z.fR5 + z.fR6 + z.fR7 + z.fR8 + z.fR9 + z.fR10 + 
	  z.fW1 + z.fW2 + z.fW3 + z.fW4 + z.fW5 + z.fW6 + 
	  z.fW7 + z.fW8 + z.fW9 + z.fW10 + z.mR1 + z.mR2 +
	  z.mR3 + z.mR4 + z.mR5 + z.mR6 + z.mR7 + z.mR8 +
	  z.mR9 + z.mR10 + z.mW1 + z.mW2 + z.mW3 + z.mW4 +
	  z.mW5 + z.mW6 + z.mW7 + z.mW8 + z.mW9 + z.mW10 - 1
	# THIS TAKES AROUND 10-15 MINUTES TO RUN
	bhbm.samples <- bhbm( data = salBhbm, 
	  fixed.formula = fixedEff, 
	  random.formula = randomEff,
	  prior = salPrior,
	  sampler = salSampler, 
	  trials.formula = ~trials,
	  overdispersion = "logit-normal" )
	
	
	# Prepare the data for a BUGS analysis
	salBugs <- list( y = salamanders$y, 
	  x = salamanders$x )
	salBugs$y <- salBugs$y[,1]
	salBugs$N <- length(salBugs$y)
	# create an indicator for the female
	salBugs$fem <- as.vector( salamanders$z[,1:20] %*% (1:20) )
	# create an indicator for the male
	salBugs$mal <- as.vector( salamanders$z[,21:40] %*% (1:20) )
	
	bugsModel <- function (){
	  for( i in 1 : N ) {
			y[i] ~ dbern(p[i])
			b[i] ~ dnorm(0.0,tau)
			logit(p[i]) <- x[i,1]*gam1 + x[i,2]*gam2 + x[i,3]*gam3 + x[i,4]*gam4 + betaFem[fem[i]] + betaMal[mal[i]] + b[i]
		}
		for( j in 1:20 ) {
		  betaFem[j] ~ dnorm(0, tauRand)
		  betaMal[j] ~ dnorm(0, tauRand)
		}		
		gam1 ~ dnorm(0.0,1)
		gam2 ~ dnorm(0.0,1)
		gam3 ~ dnorm(0.0,1)
		gam4 ~ dnorm(0.0,1)
		tau ~ dgamma(1.5, 1.5)
		sigma <- 1/sqrt(tau)
		tauRand ~ dgamma(1.5, 1.5)
		sigmaRand <- 1/sqrt(tauRand)
		invisible()
	}
	
	bugsInits <- list( gam1 = 0, gam2 = 0, gam3 = 0,
	  gam4 = 0, tau = 1, tauRand = 1 )
	bugsInits <- list( bugsInits )
	
	bugs.samples <- posteriorSamples( data = salBugs,
	  parametersToSave = c( "gam1", "gam2", "gam3", 
	  "gam4", "sigma", "sigmaRand", "betaFem", 
	  "betaMal", "p" ), 
	  model = bugsModel, 
	  nIter = 1000, inits = bugsInits,
	  nThin = 50, DIC = F, engine = "OpenBUGS" )
	
	# compare the sample distributions
	compareSampleDistributions(getSamples(bugs.samples), 
	  getSamples(bhbm.samples), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100)) 

}
