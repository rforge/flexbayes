{cat("----- No extraneous info printed while fitting the salamanders data --\n");T}
{
# Functions: bhbm
# Description: Make sure that no extraneous info is printed 
# while fitting a model for the salamanders data
  
  # save printed output to a file
  sink( "temp.txt" )
  
  # Specify the prior
  salPrior <- bhbm.prior ( 
    sigma2 = bayes.invChisq(df = 2, sigma0.sq = 10),
    fixed.coef = bayes.normal( mean = rep(0, 4), 
    cov = rep(10, 4) ),
    random.var = bayes.invChisq(df = 2, sigma0.sq=10),
    common.glm = 2 )

  # Specify the sampling parameters  
  salSampler <- bhbm.sampler( nSamples = 1000,
    nThin = 10 )

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
  # THIS TAKES A FEW MINUTES TO RUN
  salPost <- bhbm( data = salBhbm, 
    fixed.formula = fixedEff, 
    random.formula = randomEff,
    prior = salPrior,
    sampler = salSampler, 
    trials.formula = ~trials,
    overdispersion = "logit-normal" )

  # perform convergence diagnostics.  Note that by 
  # default only the first thirty parameters are 
  # displayed.
  e <- effectiveSize(salPost)

  # summarize the inferences.  Again only the first
  # thirty parameters are displayed.  
  s <- summary(salPost)

  sink()
  ( file.lengths("temp.txt")[1] == 0 )
}
