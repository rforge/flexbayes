\name{FlexBayes}
\alias{FlexBayes}
\title{Bayesian model-fitting}

\description{
  Examples of fitting Bayesian models using FlexBayes are available in the following 
	help files:
	\link{epilepsy}, \link{crossoverTrial}, \link{drugAdverseEvents}, 
	\link{effToxDose}, \link{ratsWeight}, \link{seeds}, \link{schools}, 
	\link{pumps}, \link{salamanders}, and \link{stackBayes}.
	The first four of these are clinical examples.  The others give analyses of well-known 
	data sets, for comparative purposes.  
	Additional examples are available in the help files for the model-fitting functions 
	\link{bhlm}, \link{bhbm}, \link{bhpm}, and \link{posteriorSamples}.

  Validation of the model-fitting functions in FlexBayes can be performed via a call to the function \code{validate.package}.  
  This is described in the user manual; the call is:
  \code{validate.package( "FlexBayes" )}

  Certain model-fitting and convergence diagnosis functions in FlexBayes require the installation of the 
  open-source Spotfire S+ packages BRugs, R2WinBUGS, and/or coda and the open-source standalone software WinBUGS.  
  If you call a function that requires one of these items, and it is not installed, you will see an error
  message informing you of the required dependencies.
  
  The Spotfire S+ packages BRugs, R2WinBUGS, and coda can be obtained from TIBCO's CSAN website at:
  http://spotfire.tibco.com/csan
  
  Copy the Windows binary versions of the necessary package(s) into one of the Spotfire S+ library directories.
  The open-source standalone software WinBUGS can be obtained from 
  http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/contents.shtml
  
  Version 1.4.3 or later of WinBUGS is required in order to use the \code{posteriorSamples} function with the 
  \code{engine = "WinBUGS"} option.  Install the latest version of WinBUGS, as well as any patches, then 
  complete free registration of the WinBUGS software; full instructions
  are available on the WinBUGS website given above.
}

\examples{
# The data frame stack.dat combines the stack.loss 
# and stack.x data sets included with Spotfire S+.  These 
# data are from the operation of a plant for the 
# oxidation of ammonia to nitric acid, measured on 21 
# consecutive days. 
# The goal is to model Loss (percent of ammonia lost 
# times 10), as a linear function of Air.Flow (air 
# flow to the plant), Water.Temp (cooling water inlet 
# temperature), and Acid.Conc (acid concentration as 
# a percentage). 

# Specify the hyperparameters
beta.mean <- 0
beta.var <- 1e6
# the following hyperparameter values result in a 
# very diffuse prior for the residual variance, one
# that has infinite mean and variance.
sigma2df <- 1  
sigma2scale <- 1

# We will fit the model using the blm function.

# Specify the prior distributions
coefPrior <- bayes.normal(
  mean.vector=rep(beta.mean,4), 
  covmat = diag(rep(beta.var, 4)))
  
varPrior <- bayes.invChisq(df=sigma2df, 
  sigma0.sq=(sigma2scale^2))
  
stackPrior <- bhlm.prior( fixed.coef = 
  coefPrior, error.var = varPrior)
  
# specify the sampling method and operating 
# characteristics
stackSampler <- bhlm.sampler( nChains=3, 
  init.point = "prior", nSamples=10000 )

# Obtain the posterior samples.  Regress
# Loss on all of the other variables in 
# the data set.
stackPost <- bhlm( fixed.formula = Loss ~ ., 
  data = stack.dat, prior = stackPrior,
  sampler = stackSampler)

# Perform convergence diagnostics
# Check the autocorrelation and effective
# sample sizes
autocorr.plot(stackPost)
effectiveSize(stackPost)

# create traceplots of the parameters.  Thin the 
# samples first to give clearer plots
thinned.samples <- window(stackPost, 
  thin = 5)
traceplot(thinned.samples)

# Perform the Geweke, Gelman-Rubin, and Heidelberger-
# Welch diagnostics.  See the help files for these 
# functions for the interpretation; here, all of the 
# diagnostics give satisfactory results.
geweke.diag(stackPost)
gelman.diag(stackPost)

# since the convergence diagnostics look good, 
# summarize the posterior distributions of the 
# parameters.  Obtain posterior means, standard 
# deviations, credible intervals, and posterior 
# density estimates.  The coefficient estimates
# are very close to the least squares estimates
# for the same model.
summary(stackPost)
densplot(stackPost)

# Obtain a "Bayesian p-value" for the coefficient for
# Air.Flow (here this means the posterior probability 
# that this coefficient is negative).  The Bayesian
# p-value is very small, so the data indicates that 
# this coefficient is greater than 0.
betaSamples <- 
  getSamples(stackPost,"Air.Flow")
sum(betaSamples < 0) / length(betaSamples)

}
