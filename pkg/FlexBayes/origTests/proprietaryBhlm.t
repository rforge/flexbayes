{cat("----- Make sure that bhlm does not load open-source libraries ------\n");T}
{
# Functions: bhlm
# Description: Make sure that no open-source libraries are loaded when calling bhlm

  if( "coda" %in% search() )
    detach("coda")
  if( "R2WinBUGS" %in% search() ) 
    detach("R2WinBUGS")
  if( "BRugs" %in% search() ) 
    detach("BRugs")

  ##### Fit a model to the Orthodont data
  # that comes with Spotfire S+.  
  ortho <- as.data.frame( 
    Orthodont[ Orthodont$Sex == "Female", ] )

  #the likelihood
  ortho.lkhd <- bhlm.likelihood( type = "normal" )

  #the priors
  error.var.nu <- 3
  error.var.sigma02 <- 0.25
  error.var <- bayes.invChisq( df = error.var.nu, 
    sigma0.sq = error.var.sigma02 ) 

  alpha.mean <- zero
  alpha.Cov <- identity
  level2.coef <- bayes.normal( 
    mean.vector = alpha.mean, covmat = alpha.Cov )

  random.var.nu <- 3
  random.var.sigma02 <- 1
  random.var <- bayes.invChisq( df = random.var.nu, 
    sigma0.sq = random.var.sigma02 ) 

  ortho.prior <- bhlm.prior( 
    error.var = error.var, random.var = random.var,
    level2.coef = level2.coef )

  #the sampler parameters
  ortho.sampler <- bhlm.sampler( init.point = "prior" )

  #the call to bhlm
  ortho.bhlm <- bhlm( 
    random.formula = distance ~ I(age - 11), 
    level2.formula = ~ 1,
    group.formula = ~ Subject, data = ortho, 
    prior = ortho.prior, 
    likelihood = ortho.lkhd, sampler = ortho.sampler )

  ( !( "coda" %in% search() ) && !( "R2WinBUGS" %in% search() ) &&
    !( "BRugs" %in% search() ) )
}
