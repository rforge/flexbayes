{cat("----- No extraneous info printed while fitting the schools model --\n");T}
{
# Functions: bhlm, posteriorSamples
# Description: Make sure that no extraneous info is printed 
# while fitting a model for the schools data
  
  # save printed output to a file
  sink( "temp.txt" )
  
  # The schools data set is found in the library 
  # R2WinBUGS
  library(R2WinBUGS)
  data(schools)

  ########################################
  # First fit a linear model using bhlm.
  # Specify the prior:
  nData <- length(schools$sd)
  error.var <- vector( "list", nData )
  for (i in (1:nData)){
    error.var[[i]] <- bayes.massPoint( schools$sd[i]^2 )
  }
  schoolsPrior <- bhlm.prior (
    error.var = error.var,
    fixed.coef = "non-informative",
    common.error.var = 0 )

  schoolsSampler <- bhlm.sampler ( nSamples = 10000,
    nThin = 10 )

  # Fit the model
  schoolsPost <- bhlm( fixed.formula = estimate~1,
    data = schools, group.formula = ~school, 
    prior = schoolsPrior,
    sampler = schoolsSampler )

  # summarize.  The results are bogus here; something is 
  # wrong with the code to fix the outcome variance, 
  # needs to be debugged.
  s <- summary( schoolsPost )

  #############################################
  # For slightly more flexible prior specification,
  # fit using the posteriorSamples function.

  # First, specify the model
  bugsModel <- function() {
    for (j in 1:J){
      y[j] ~ dnorm (theta[j], y.prec[j])
      y.prec[j] <- pow(y.sigma[j], -2)
      theta[j] ~ dnorm (theta.mean, theta.prec)
    }
    theta.mean ~ dunif (-100, 100)
    theta.prec <- pow (theta.sigma, -2)
    theta.sigma ~ dunif (0, 100)
  }  

  # Specify the data set 
  data <- list( y = schools$estimate, 
          y.sigma = schools$sd, 
          J = length(schools$estimate) )

  # Specify the initial values
  inits <- list(theta.sigma = 50, theta.mean = 0)

  # There need to be initial values specified for each
  # chain that will be run, so the inits object needs 
  # to be a list of initial value lists.  Here we will 
  # only run one chain.
  inits <- list(inits)

  # Specify which parameters should be stored during the
  # model run for use in inference.
  paramsToSave <- c("theta", "theta.mean", "theta.sigma")

  # Run the model and obtain the posterior samples
  schoolsPost <- posteriorSamples(model = bugsModel, 
    data = data, inits = inits, 
    parametersToSave = paramsToSave, 
    nIter = 10000, engine = "WinBUGS")

  a <- autocorr(schoolsPost)
  c <- crosscorr(schoolsPost)
  e <- effectiveSize(schoolsPost)
  g <- geweke.diag(schoolsPost)
  h <- heidel.diag(schoolsPost)

  # Run the chain for more iterations and use thinning
  # in order to decrease the autocorrelation of the 
  # parameters.  Keep the burnin length short because 
  # there appeared from the traceplots to have been enough
  # burnin.  Another option to improve mixing could be to 
  # reducing the crosscorrelation by reparameterization of 
  # the model.  

  schoolsPost <- posteriorSamples(model = bugsModel, 
    data = data, inits = inits, 
    parametersToSave = paramsToSave,
    nIter = 10000, nThin = 100, engine = "WinBUGS")

  e <- effectiveSize(schoolsPost)
  # The above diagnostics all indicate good convergence 
  # and mixing of the chain.  The effective sample sizes 
  # are all at least 5000.

  c <- crosscorr(schoolsPost)
  s <- summary(schoolsPost)

  sink()
  ( file.lengths("temp.txt")[1] == 0 )
}
