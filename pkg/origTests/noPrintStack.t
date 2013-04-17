{cat("----- No extraneous info printed while fitting the stacks data --\n");T}
{
# Functions: bhlm, posteriorSamples
# Description: Make sure that no extraneous info is printed 
# while fitting a model for the stacks data
  
  # save printed output to a file
  sink( "temp.txt" )
  
  # We fit a linear model of Loss with Air.Flow, 
  # Water.Temp, and Acid.Conc as covariates

  # Specify the hyperparameters
  beta.mean <- 0
  beta.prec <- 1e-6
  # the following hyperparameter values result in a 
  # very diffuse prior for the residual variance, one
  # that has infinite mean and variance.
  sigma2df <- 1  
  sigma2scale <- 1

  # First we fit the model using the bhlm function.

  # Specify the prior distributions
  coefPrior <- bayes.normal(
    mean.vector=rep(beta.mean,4), 
    covmat = diag(rep(beta.prec^-1, 4)))
  
  varPrior <- bayes.invChisq(df=sigma2df, 
    sigma0.sq=(sigma2scale^2))
  
  stackPrior <- bhlm.prior( fixed.coef = coefPrior, 
    error.var = varPrior)
  
  # specify the sampling method and other 
  # operating characteristics for the sampler
  stackSampler <- bhlm.sampler( nChains=3, 
    init.point = "prior", 
    nSamples=10000 )

  # Obtain the posterior samples.  Regress
  # Loss on all of the other variables in 
  # the data set.
  stackSamples <- bhlm( fixed.formula = Loss ~ ., 
    data = stack.dat, 
    prior = stackPrior,
    sampler = stackSampler)

  # Perform convergence diagnostics
  g <- geweke.diag(stackSamples)
  g <- gelman.diag(stackSamples)

  # Since the convergence diagnostics are 
  # satisfactory, summarize the results.
  s <- summary(stackSamples)

  # As an alternative, this model can be fit
  # using the WinBUGS engine.  This allows for
  # greater flexibility but requires a more 
  # detailed specification.
  ## specify the model in BUGS format     
  bugsModel <- function (){
    for( i in 1 : N ) {
      Loss[i] ~ dnorm(mu[i],tau.c)
      mu[i] <- beta[1] + beta[2] * Air.Flow[i] + 
        beta[3]*Water.Temp[i] + beta[4]*Acid.Conc.[i]
    }
    for(j in 1:4){
      beta[j] ~ dnorm(beta.mean, beta.prec)
    }
    tau.c ~ dgamma(a.prec, b.prec)
    sigma <- pow(tau.c, -0.5)

    invisible()           
  }

  # Create the data frame for fitting the model
  stackModelData <- stack.dat
  stackModelData$N <- length(stack.dat$Loss)

  # Set the values for the hyperparameters.  Choose 
  # values that lead to fairly diffuse priors on the 
  # variance 1/tau.c and the coefficients beta[j]
  a.prec <- 0.5
  b.prec <- 1
  stackModelData$a.prec <- a.prec
  stackModelData$b.prec <- b.prec
  stackModelData$beta.mean <- beta.mean
  stackModelData$beta.prec <- beta.prec

  # Specify the parameters for which we wish to save 
  # the posterior samples
  parameters.to.save <- c("beta", "sigma")

  # Run multiple chains, so that the Gelman-Rubin 
  # diagnostic can be used.
  nChains <- 3

  # Make a list of initial values for each of the 
  # chains.  Draw the initial values from the prior 
  # distribution.  Since the prior distributions are 
  # fairly vague, these values should be over-dispersed 
  # relative to the posterior distribution.  This will 
  # allow us to use the Gelman-Rubin convergence 
  # diagnostic.
  initialValues <- list(rep(-1, nChains))
  for (i in (1:nChains)){
    betaInit <- rnorm(n = 4, mean = beta.mean, 
      sd = 1 / sqrt(beta.prec))
    tauInit <- rgamma(n = 1, shape = a.prec, 
      rate = b.prec) 
    initsThisChain <- list(beta = betaInit, 
      tau.c = tauInit)
    initialValues[[i]] <- initsThisChain 
  }

  ## obtain the posterior samples
  posterior.samples <- posteriorSamples (
    data = stackModelData, model = bugsModel, 
    inits = initialValues, nChains = nChains,
    parametersToSave = parameters.to.save,
    nIter = 10000, nBurnin = 1000, nThin = 5,
    engine = "WinBUGS")

  ## perform convergence diagnostics
  e <- effectiveSize(posterior.samples)

  # perform the Gelman-Rubin diagnostic.  The 
  # 'potential scale reduction factors' should be close 
  # to 1.  Gelman et al. suggest that any value less 
  # than 1.2 is acceptable for most examples.  Here it 
  # is very close to 1.  Transform the parameters, 
  # since the Gelman-Rubin diagnostic is based on a 
  # normal approximation to the posterior distribution, 
  # and one of our parameters is a variance parameter.
  g <- gelman.diag(posterior.samples, transform = T)

  # perform the Heidelberger-Welch convergence tests, 
  # which pass.
  h <- heidel.diag(posterior.samples)

  # obtain the Geweke Z-scores.  They are mostly 
  # between -2 and 2, so there is no evidence of lack 
  # of convergence here.
  g <- geweke.diag(posterior.samples)

  # since the convergence diagnostics look good, 
  # summarize the posterior distributions of the 
  # parameters.  Obtain posterior means, standard 
  # deviations, credible intervals, and posterior 
  # density estimates.
  s <- summary(posterior.samples)

  # Obtain a "Bayesian p-value" for beta[2] (the 
  # posterior probability that beta[2] < 0).  It is 
  # very small, so the data indicates that beta[2] is 
  # greater than 0.
  betaSamples <- 
    getSamples(posterior.samples,"beta[2]")
  s <- sum(betaSamples < 0) / length(betaSamples)
  
  sink()
  ( file.lengths("temp.txt")[1] == 0 )
}
