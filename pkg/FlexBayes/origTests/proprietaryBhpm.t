{cat("----- Make sure that bhpm does not load open-source libraries ------\n");T}
{
# Functions: bhpm
# Description: Make sure that no open-source libraries are loaded when calling bhpm

  if( "coda" %in% search() )
    detach("coda")
  if( "R2WinBUGS" %in% search() ) 
    detach("R2WinBUGS")
  if( "BRugs" %in% search() ) 
    detach("BRugs")

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

  ( !( "coda" %in% search() ) && !( "R2WinBUGS" %in% search() ) &&
    !( "BRugs" %in% search() ) )
}
