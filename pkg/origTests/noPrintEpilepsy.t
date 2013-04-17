{cat("----- No extraneous info printed while fitting epilepsy model --\n");T}
{
# Functions: bhpm
# Description: Make sure that no extraneous info is printed 
# while fitting a model to the epilepsy data
  
  # save printed output to a file
  sink( "temp.txt" )
  
  # Specify the prior
  epilepsyPrior <- bhpm.prior ( 
    sigma2 = bayes.nonInfoPower(-1), 
    common.glm = 2)

  # Specify the sampling parameters  
  epilepsySampler <- bhpm.sampler( nSamples = 1000,
    nThin = 10, init.point = "user's choice",
    params.init = list( sigma2 = 1, 
    fixed.coef = rep(0, 1), random.var = 1,
    random.coef = 0, level2.coef = rep(0, 5) ) )

  # Load the data
  data(epilepsy)
  # Prepare the data 
  nSubj <- dim(epilepsy$y)[1]
  base <- log ( epilepsy$baseline / 4 )
  epilepsyBhpm <- data.frame( 
    y = as.vector( epilepsy$y ),
    subj = rep( (1:nSubj), each = 4 ),
    visit4 = rep( epilepsy$visit4, nSubj ),
    base = rep( base, each = 4 ),
    treat = rep( epilepsy$treatment, each = 4 ),
    logAge = rep( log( epilepsy$age ), each = 4 ) )  

  # fit the model
  fixedEff <- y ~ visit4 - 1
  level2Eff <- ~ treat + base + logAge + treat * base
  randomEff <- ~ 1
  # THIS TAKES A FEW MINUTES TO RUN
  epilepsyPost <- bhpm( data = epilepsyBhpm, 
    fixed.formula = fixedEff, 
    random.formula = randomEff,
    group.formula = ~ subj,
    level2.formula = level2Eff,
    prior = epilepsyPrior,
    sampler = epilepsySampler, 
    overdispersion = "log-normal" )

  g <- geweke.diag(epilepsyPost)

  # The convergence diagnostics look good, so 
  # summarize the inferences.  Again only the first
  # thirty parameters are displayed.  The (Intercept):#
  # parameters are the random effects.  
  # The results are similar to those reported by 
  # Breslow and Clayton (1993).
  s <- summary(epilepsyPost)
  # Look at the names of the rest of the parameters in 
  # the model.  The lambda parameters are data 
  # augmentation parameters for the overdispersion
  # (one lambda parameter for each data point).
  v <- varnames(epilepsyPost)
  
  sink()
  ( file.lengths("temp.txt")[1] == 0 )
}
