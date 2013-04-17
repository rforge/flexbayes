{cat("----- Make sure that bhpm does not print extraneous information ------\n");T}
{
# Functions: bhpm
# Description: Make sure that bhpm does not print extraneous information
  
  # save printed output to a file
  sink( "temp.txt" )
  
  ##########################################
  # Use bhpm to fit an overdispersed Poisson model to the
  # pumps data

  xi.prior = bayes.uniformShrinkage (0.5)
  random.var.prior = bayes.invWishart( df = 3, scale = diag( c(1,1) ) )

  pump.prior = bhpm.prior ( xi = xi.prior,
    random.var =  random.var.prior, 
    common.glm = 2 )

  # Specify the control parameters (burn-in length, 
  # etc.) and the initial values for the parameters.  
  pump.sampler <- bhpm.sampler( nSamples=1000, nThin= 50, 
    nChains = 1, nBurnin = 1000,
    init.point = "prior", update.cov = 0 ) 

  pump.exposure <- ~ e
  pump.random <- z ~ 1 + x

  ## call bhpm to fit a gamma-conjugate model
  pump.bhpm <- bhpm( random.formula = pump.random, 
    exposure.formula = pump.exposure, data = pumps, 
    prior = pump.prior, sampler = pump.sampler,
    overdispersion = "gamma-conj" )

  # summarize the posterior distribution
  s <- summary(pump.bhpm)

  ##########################################
  # Use bhpm to fit a hierarchical model to
  # the epilepsy data

  # Specify the prior
  epilepsyPrior <- bhpm.prior ( 
    sigma2 = bayes.invChisq(df = 3, sigma0.sq = 10),
    fixed.coef = bayes.normal( mean = 0, cov = 10 ),
    level2.coef = bayes.normal( mean = rep(0, 5), cov = diag( rep(10, 5) ) ),
    random.var = bayes.invChisq(df = 3, sigma0.sq = 10),
    common.glm = 2 )

  # Specify the sampling parameters  
  epilepsySampler <- bhpm.sampler( nSamples = 1000,
    nThin = 10 )

  # Load the data
  data(epilepsy)
  # Prepare the data 
  nSubj <- dim(epilepsy$y)[1]
  base <- log ( epilepsy$baseline / 4 )
  epilepsyBhpm <- data.frame( 
    y = as.vector( t( epilepsy$y ) ),
    subj = rep( (1:nSubj), each = 4 ),
    visit4 = rep( epilepsy$visit4, nSubj ),
    base = rep( base, each = 4 ),
    treat = rep( epilepsy$treatment, each = 4 ),
    logAge = rep( log( epilepsy$age ), each = 4 ) )  

  # fit the model
  fixedEff <- 
    y ~ visit4 - 1
  level2Eff <- ~ treat + base + logAge + treat * base
  randomEff <- ~ 1
  # THIS TAKES A FEW MINUTES
  epilepsyPost <- bhpm( data = epilepsyBhpm, 
    fixed.formula = fixedEff, 
    random.formula = randomEff,
    group.formula = ~ subj,
    level2.formula = level2Eff,
    prior = epilepsyPrior,
    sampler = epilepsySampler, 
    overdispersion = "log-normal" )

  # The convergence diagnostics look good, so 
  # summarize the inferences
  s <- summary(epilepsyPost)
  
  sink()
  ( file.lengths("temp.txt")[1] == 0 )
}
