{cat("----- Make sure that bhbm does not load open-source libraries ------\n");T}
{
# Functions: bhbm
# Description: Make sure that no open-source libraries are loaded when calling bhbm

  if( "coda" %in% search() )
    detach("coda")
  if( "R2WinBUGS" %in% search() ) 
    detach("R2WinBUGS")
  if( "BRugs" %in% search() ) 
    detach("BRugs")

  # Fit an overdispersed binomial model to the 
  # toxoplasmosis data (Kass and Steffey 1989).
  toxo.prior <- bhbm.prior( xi = 
    bayes.uniformShrinkage( 1 ), 
    random.var = bayes.invChisq( df=3, sigma0.sq=1 ),
    common.glm = 2 )

  # Specify the control parameters (burn-in length, 
  # etc.) and the initial values for the parameters.  
  # The specification of an initial value for a
  # fixed effect is required, but the value will be 
  # ignored because no fixed effects are specified in 
  # the call to bhbm.
  toxo.sampler <- bhbm.sampler( nBurnin = 10000, 
    nSamples = 1000, nThin = 50, nChains = 1, 
    update.cov = 1, init.point = "user's choice", 
    params.init = list( xi = 10, random.coef = 0.5, 
    fixed.coef = 0, random.var = 1 ) )  

  toxo.trials <- ~ ni
  toxo.random <- yi ~ 1

  ##########################################
  ## call bhbm to fit a beta-conjugate model
  toxo.bhbm <- bhbm( random.formula = toxo.random, 
    trials.formula = toxo.trials, data = toxo.dat, 
    overdispersion = "beta-conj", 
    prior = toxo.prior, sampler = toxo.sampler )

  ( !( "coda" %in% search() ) && !( "R2WinBUGS" %in% search() ) &&
    !( "BRugs" %in% search() ) )
}
