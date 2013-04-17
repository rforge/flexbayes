{cat("----- No extraneous info printed while fitting the pumps data --\n");T}
{
# Functions: bhpm, posteriorSamples
# Description: Make sure that no extraneous info is printed 
# while fitting a model for the schools data
  
  # save printed output to a file
  sink( "temp.txt" )
  
  # First fit the pumps data using the bhpm function

  # Specify the priors
  pump.prior = bhpm.prior ( xi = 
    bayes.uniformShrinkage (0.225), 
    fixed.coef = "non-informative",
    common.glm = 2 )

  # Specify the chain controls (burn-in length, 
  # etc.) and the initial values for the parameters.
  pump.sampler <- bhpm.sampler( nSamples=1000, 
    nThin= 50, nChains = 3, nBurnin = 1000,
    init.point = "prior", update.cov = 1 ) 

  ##########################################
  ## call bhpm to fit a gamma-conjugate model
  pump.exposure <- ~ e
  pump.fixed <- z ~ x
  pump.bhpm <- bhpm( fixed.formula = pump.fixed, 
    exposure.formula = pump.exposure, data = pumps, 
    prior = pump.prior, sampler = pump.sampler,
    overdispersion = "gamma-conj" )

  # check convergence
  g <- geweke.diag(pump.bhpm)

  # Summarize the posterior distribution.
  # The results are close to those reported 
  # by Christiansen and Morris (1997).
  s <- summary(pump.bhpm)

  ##############################################
  # Fit a similar model on the pumps data using 
  # the posteriorSamples function.

  # First specify a Poisson model using BUGS
  # syntax.
  # The outcome z[i] for each pump is the number of 
  # failures of the pump.  The number of failures is 
  # modeled via a Poisson distribution with mean equal 
  # to e[i]*lambda[i], where e[i] is the duration of 
  # operation in units of 1,048 hours and lambda[i] is 
  # the failure rate of the pump per unit time.  The 
  # parameters lambda[i] are modeled as random effects 
  # by giving them a common prior distribution for all 
  # of the pumps.  A conjugate gamma distribution is 
  # chosen, with mean mu[i] and variance mu[i]^2/xi,
  # so that xi is a precision parameter.  mu[i] is then
  # regressed on a predictor variable, which is an 
  # indicator of whether the pump was operated 
  # continuously (x[i] = 1) or intermittently 
  # (x[i] = 0).

  # First we illustrate a mistaken model specification,
  # and how to diagnose the problem.  

  bugsModel <- function() {
    for (i in 1:N){
      z[i] ~ dpois(gamma[i])
      gamma[i] <- e[i] * lambda[i]
      lambda[i] ~ dgamma(xi, b[i])
      b[i] <- xi/mu[i]
      mu[i] <- beta0 + beta1*x[i]
        # The mistake here is to set mu[i] equal to the
        # linear regression function instead of 
        # using a log link function as in the 
        # original model specification.
     }
    xi <- z0*(1-B0)/B0
      # There is a second mistake, namely the use of
      # the constant z0, the value of which has not 
      # been specified.
    B0 ~ dunif(0,1)
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    invisible()
  }

  pumpsBugs <- list( z = pumps$z, e = pumps$e, 
    x = pumps$x, N = length(pumps$z) )

  inits <- list(beta0 = 0, beta1 = 0, B0 = 0.5)
  inits <- list(inits)

  # Add the value of z0 to the data set.  Here we 
  # choose the value used by Christiansen and Morris 
  # 1997.
  pumpsBugs <- c(pumpsBugs, z0=0.225)

  # Inspection of the model and the initial values 
  # reveals that there is a division by zero for the 
  # initial values that we chose, due to the model
  # misspecification.  Fix the model.
  bugsModel <- function() {
    for (i in 1:N){
      z[i] ~ dpois(gamma[i])
      gamma[i] <- e[i] * lambda[i]
      lambda[i] ~ dgamma(xi, b[i])
      b[i] <- xi/mu[i]
      mu[i] <- exp(beta0 + beta1*x[i])
    }
    xi <- z0*(1-B0)/B0
    B0 ~ dunif(0,1)
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    invisible()
  }

  # Now that the model is running correctly, run it
  # in WinBUGS within Spotfire S+.
  pumpsPost <- posteriorSamples(model = bugsModel, 
    data = pumpsBugs, inits = inits, 
    parametersToSave = 
      c("beta0", "beta1", "xi", "lambda"),
    nIter = 10000, nThin = 10, engine = "WinBUGS")

  e <- effectiveSize(pumpsPost)
  # Perform the Geweke diagnostic.  The resulting values 
  # are z-scores, which in this case are not extreme.
  g <- geweke.diag(pumpsPost)

  # Perform the Heidelberger-Welch diagnostic, which 
  # passes.
  h <- heidel.diag(pumpsPost)

  # The convergence diagnostics look good, so summarize 
  # the inferences.  The lambda[i] parameters are the 
  # pump-specific failure rates (per 1,048 hours).  
  # There is a great deal of variability in these 
  # rates; the interval estimate of one of these rates
  # is far above the others.  
  # Notice that the estimates of the regression 
  # coefficients here are close to those in 
  # Christiansen and Morris (1997), Table 1.,
  # and identical to the results obtained from the
  # bhpm function above.
  s <- summary(pumpsPost)

  sink()
  ( file.lengths("temp.txt")[1] == 0 )
}
