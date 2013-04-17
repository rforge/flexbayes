\name{schools}
\alias{schools}
\docType{data}
\title{8 schools analysis}
\description{The 8 schools data set from Rubin (1981).  In the example
  below we run the model given in Gelman, Carlin, Stern,
  and Rubin (2004), 
  perform convergence analysis, and summarize the posterior inferences.}
\usage{data(schools)}
\format{
  A data frame with 3 variables and 8 observations. The variables are:
  \describe{
    \item{school}{See Source.}
    \item{estimate}{See Source.}
    \item{sd}{See Source.}
  }
}
\source{
Rubin, D.B. (1981): Estimation in Parallel Randomized Experiments.
\emph{Journal of Educational Statistics} 6(4), 377-400.

Section 5.5 of
Gelman, A., Carlin, J.B., Stern, H.S., Rubin, D.B. (2004): 
\emph{Bayesian Data Analysis}, 2nd edition, CRC Press.
}
\keyword{datasets}

\examples{

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

# check convergence
autocorr.plot( schoolsPost )

# summarize.  The results are bogus here; something is 
# wrong with the code to fix the outcome variance, 
# needs to be debugged.
summary( schoolsPost )

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

traceplot(schoolsPost)
# The traceplot reveals some poor mixing of the 
# chain.

autocorr.plot(schoolsPost)
# The autocorrelation plots show some autocorrelation 
# of the chain up through lag 30.

autocorr(schoolsPost)
# The lag-50 autocorrelation of mu.theta and of some 
# of the theta parameters is over 5 percent

crosscorr(schoolsPost)
# There is close to 50 percent cross-correlation between
# some pairs of parameters.  
# This is probably causing the slow mixing of the chain.  

effectiveSize(schoolsPost)
# The effective sample sizes of some of the parameters 
# are less than 500, although the MCMC drew 10,000 
# samples.  This is due to the autocorrelation.

geweke.diag(schoolsPost)
geweke.plot(schoolsPost)
# The Geweke diagnostic finds the problem, as do  
# the Geweke-Brooks plots.  For the plots, many
# of the Z-scores are outside of the 95-percent 
# intervals.

heidel.diag(schoolsPost)
# The halfwidth tests in the Heidelberger-Welch 
# diagnostic fail for some of the parameters.  The 
# stationarity tests all pass, despite the 
# autocorrelation problem.

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

# create traceplots of the parameters.  Increase the 
# thinning of the samples first to give clearer plots.
thinned.samples <- window(schoolsPost, thin = 1000)
traceplot(thinned.samples)

autocorr.plot(schoolsPost)
effectiveSize(schoolsPost)
# The above diagnostics all indicate good convergence 
# and mixing of the chain.  The effective sample sizes 
# are all at least 5000.

geweke.plot(schoolsPost)
# The Geweke plots are acceptable.  There are a few 
# Z-scores that fall outside of the 95-percent 
# intervals; however, we expect this due to multiple 
# testing.  As long as there are not many falling 
# outside of the intervals, the plots are considered 
# acceptable.

crosscorr(schoolsPost)
# Note that the cross-correlations of the parameters 
# are still high.  Increasing the thinning of the chain 
# does not affect the cross-correlations.  However, 
# having high cross-correlation is not a problem if the 
# autocorrelations are low.

# Since the convergence diagnostics look good, summarize 
# the posterior inferences.  The results here match those
# given in Gelman, Carlin, Stern, and Rubin (2004).
summary(schoolsPost)

}
