\name{pumps}
\alias{pumps}
\docType{data}
\title{Failure data for power plant pumps}
\description{Failure data for power plant pumps, as analyzed in Christiansen
and Morris (1997).  Here we fit a binomial model, assess convergence, and obtain posterior inferences.}
\usage{data(pumps)}
\source{
Christiansen, CL and Morris, CN (1997).  Hierarchical Poisson Regression Modeling.
Journal of the American Statistical Association, 92, 618-632.
}
\keyword{datasets}

\examples{

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
autocorr.plot(pump.bhpm)
traceplot(pump.bhpm)
geweke.diag(pump.bhpm)

# Summarize the posterior distribution.
# The results are close to those reported 
# by Christiansen and Morris (1997).
pump.bhpm

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

# The first time a model is run, we suggest using
# the WinBUGS engine along with the debug=T option,
# which allows for viewing of the BUGS error 
# messages.  After viewing the error messages,
# close the WinBUGS window manually.
# (Ignore the following "Not run" statements; these
# indicate to Spotfire S+ not to run the following
# statements while checking the package)
\dontrun{pumpsPost <- posteriorSamples(model = 
  bugsModel, data = pumpsBugs, inits = inits, 
  parametersToSave = 
    c("beta0", "beta1", "xi", "lambda"),
  nIter = 10000, nThin = 10, engine = "WinBUGS",
  debug = T)
}
# Add the value of z0 to the data set.  Here we 
# choose the value used by Christiansen and Morris 
# 1997.
pumpsBugs <- c(pumpsBugs, z0=0.225)

# Run the model again, which shows that there is 
# still a numerical problem.
\dontrun{pumpsPost <- posteriorSamples(model = bugsModel, 
  data = pumpsBugs, inits = inits, 
  parametersToSave = 
    c("beta0", "beta1", "xi", "lambda"),
  nIter = 10000, nThin = 10, engine = "WinBUGS",
  debug = T)
}
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

# Rerun the model, still in the debug mode.  The
# model will run correctly, but afterwards the 
# WinBUGS window will stay open.  Close it manually.
\dontrun{pumpsPost <- posteriorSamples(model = bugsModel, 
  data = pumpsBugs, inits = inits, 
  parametersToSave = 
    c("beta0", "beta1", "xi", "lambda"),
  nIter = 10000, nThin = 10, engine = "WinBUGS",
  debug = T)
}
# Now that the model is running correctly, run it
# using WinBUGS
pumpsPost <- posteriorSamples(model = bugsModel, 
  data = pumpsBugs, inits = inits, 
  parametersToSave = 
    c("beta0", "beta1", "xi", "lambda"),
  nIter = 10000, nThin = 10, engine = "WinBUGS")

# Apply convergence diagnostics
thinned.samples <- window(pumpsPost, thin = 100)
traceplot(thinned.samples)
effectiveSize(pumpsPost)
autocorr.plot(pumpsPost)
# Apply the Geweke diagnostic.  The resulting values 
# are z-scores, which in this case are not extreme.
geweke.diag(pumpsPost)

# Apply the Heidelberger-Welch diagnostic, which 
# passes.
heidel.diag(pumpsPost)

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
summary(pumpsPost)

}
