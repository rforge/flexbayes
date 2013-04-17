{cat("----- Make sure that bhlm does not print extraneous information ------\n");T}
{
# Functions: bhlm
# Description: Make sure that bhlm does not print extraneous information
  
  sink( "temp.txt" )
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

  s <- summary(ortho.bhlm)

  #### Fit a model without second-level coefficients
  ortho.bhlm <- bhlm(random.formula = distance ~ I(age-11), 
    fixed.formula = ~ I(age-11) + Sex + I(age-11)*Sex,
    group.formula = ~ Subject, data = Orthodont )

  s <- summary(ortho.bhlm)
  
  sink()
  ( file.lengths("temp.txt")[1] == 0 )
}
