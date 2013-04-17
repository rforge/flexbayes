{cat("----- Calls to getSamples do not return an mcmc object --\n");T}
{
# Functions: getSamples
# Description: Make sure that calls to the getSamples
# function do not return an mcmc or posterior 
# object, that instead they return a matrix or vector
# object as appropriate.  
# Use the ratsWeight example for this. 

  # Choose the hyperparameter values.  They
  # are taken to be equal to
  # those given in Gelfand et al. (1990).
  R <- diag( c( 100, 0.1 ) )
  rho <- 2

  # First, fit the model using bhlm.
  # Specify the priors.
  error.var.prior <- bayes.invChisq( df = 3, 
    sigma0.sq = 10 )
  random.var.prior <- bayes.invWishart( df = rho, 
    scale = solve( R ) )
  rats.prior <- bhlm.prior( error.var = error.var.prior,
    random.var = random.var.prior,
    level2.coef = "non-informative" )

  # Specify the control variates for the chain
  rats.sampler <- bhlm.sampler( nSamples = 1000,
    nThin = 1, init.point = "prior", nChain = 2 )

  # Prepare the data to be fit using bhlm
  rats.bhlm <- data.frame( y = as.vector( ratsWeight$y ),
    x = rep( x = ratsWeight$x, each = 30 ),
    id = rep( (1:30), 5 ) )

  # Fit the model using bhlm
  bhlm.samples <- bhlm( data = rats.bhlm,
    random.formula = y ~ x, group.formula = ~id,
    level2.formula = ~ 1,
    prior = rats.prior, sampler = rats.sampler )
  
  loadCoda()
  a <- getSamples( bhlm.samples )
  b <- getSamples( bhlm.samples, chainIndex = 2 )
  c <- getSamples( bhlm.samples, param = 2:3 )
  d <- getSamples( bhlm.samples, param = 2:3,
    chainIndex = 2 )
  e <- getSamples( bhlm.samples, param = 2 )
  f <- getSamples( bhlm.samples, param = 2,
    chainIndex = 2 )
  allMatVec <- is.matrix(a) && is.matrix(b) && is.matrix(c) &&
    is.matrix(d) && is.vector(e) && is.vector(f)
  
  noneMcmc <- !is.mcmc(a) && !is.mcmc(b) && !is.mcmc(c) &&
    !is.mcmc(d) && !is.mcmc(e) && !is.mcmc(f)
  
  correctClass <- ( class(a) == "matrix" ) && 
    ( class(b) == "matrix" ) && ( class(c) == "matrix" ) && 
    ( class(d) == "matrix" ) && ( class(e) == "numeric" ) &&
    ( class(f) == "numeric" )
  
  (allMatVec && noneMcmc && correctClass)
}
