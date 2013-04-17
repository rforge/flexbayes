{cat("----- Make sure that utility fns do not load open-source libraries --\n");T}
{
# Functions: bhlm
# Description: Make sure that no open-source libraries 
# are loaded when calling functions that manipulate 
# posterior objects, such as nvar, subsetting operations,
# etc.  Use the ratsWeight example for this.
  
  if( "coda" %in% search() )
    detach("coda")
  if( "R2WinBUGS" %in% search() ) 
    detach("R2WinBUGS")
  if( "BRugs" %in% search() ) 
    detach("BRugs")
  
  # Load the data.  This is data for the weights of 
  # infant rats during the first 36 days of life.
  # The original study compares the growth rates of
  # the rats in a control group to those in a 
  # treatment group.  Here we estimate the growth 
  # rates in the control group only.
  data(ratsWeight)

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
    nThin = 10, init.point = "prior" )

  # Prepare the data to be fit using bhlm
  rats.bhlm <- data.frame( y = as.vector( ratsWeight$y ),
    x = rep( x = ratsWeight$x, each = 30 ),
    id = rep( (1:30), 5 ) )

  # Fit the model using bhlm
  bhlm.samples <- bhlm( data = rats.bhlm,
    random.formula = y ~ x, group.formula = ~id,
    level2.formula = ~ 1, 
    prior = rats.prior, sampler = rats.sampler )
    
  a <- nvar( bhlm.samples )
  a <- varnames( bhlm.samples )
  a <- niter( bhlm.samples )
  a <- nchain( bhlm.samples )
  a <- thin( bhlm.samples )
  a <- start( bhlm.samples )
  a <- end( bhlm.samples )
  a <- summary( bhlm.samples )
  a <- 
    bhlm.samples[,c("SIGMA", "(Intercept):(Intercept)")]
  a <- bhlm.samples[[1]]
  a <- bhlm.samples[1:10,]
  a <- bhlm.samples[1:10,"SIGMA"]
  a <- as.mcmc.list( bhlm.samples )
  a <- getSamples( bhlm.samples )
  a <- getSamples( bhlm.samples, 
    param = c("SIGMA", "(Intercept):(Intercept)") )
  a <- getSamples( bhlm.samples, 
    param = c("SIGMA", "(Intercept):(Intercept)"),
    chainIndex = 1 )

  ( !( "coda" %in% search() ) && !( "R2WinBUGS" %in% search() ) &&
    !( "BRugs" %in% search() ) )
}
