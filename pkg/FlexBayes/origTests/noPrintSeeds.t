{cat("----- No extraneous info printed while fitting seed germination model --\n");T}
{
# Functions: bhbm
# Description: Make sure that no extraneous info is printed 
# while fitting a model for seed germination rates
  
  # save printed output to a file
  sink( "temp.txt" )
  
  # Load the data
  data(seeds)

  # Specify a prior that is fairly diffuse relative 
  # to this context
  seedsPrior <- bhbm.prior ( 
    sigma2 = bayes.invChisq(df = 1, sigma0.sq = 1),
    fixed.coef = bayes.normal(mean = rep(0, 4), 
    cov = diag( rep(100, 4) ) ),
    common.glm = 2 )

  # Specify the sampling parameters  
  seedsSampler <- bhbm.sampler( nSamples = 1000,
    nThin = 10 )

  # fit the model
  seedsPost <- bhbm( data = seeds, 
    trials.formula = ~n, fixed.formula = r~x1+x2+x1*x2,
    overdispersion = "logit-normal", prior = seedsPrior,
    sampler = seedsSampler )

  # perform convergence diagnostics
  g <- geweke.diag(seedsPost)

  # The convergence diagnostics look good, so 
  # summarize the inferences
  s <- summary(seedsPost)

  sink()
  ( file.lengths("temp.txt")[1] == 0 )
}
