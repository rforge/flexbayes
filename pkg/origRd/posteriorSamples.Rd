\name{posteriorSamples}
\alias{posteriorSamples}
\title{Sample from a Posterior Distribution}

\usage{
posteriorSamples(data, inits = NULL, 
  parametersToSave,  model = "model.txt", 
  nChains = 1, nIter = 10000, nBurnin = 1000,
  nThin = max(1, floor(nIter / 10000)),
  DIC = TRUE, workingDirectory = NULL, 
  bugsDirectory = getBUGSDirectory(),
  debug = FALSE, engine = "OpenBUGS", 
  overRelax = FALSE)
}
\arguments{
  \item{data}{a named list (names corresponding to variable names in the model) of the data}
  \item{inits}{a list with \code{nChains} elements; each element of the list is
    itself a list of starting values for the parameters in the model.
    Alternatively, if \code{inits = NULL}, initial values are generated}
  \item{parametersToSave}{character vector of the names of the parameters that should be monitored}
  \item{model}{either the name of the model file or the model, specified as a function.  The model should be
  	specified in BUGS format.  If it is specified in a file, the file extension can be either \file{.bug} 
  	or \file{.txt}.  For more information on specifying the model using a function, see the help for \code{write.model}
  	in the \code{R2WinBUGS} library.}
  \item{nChains}{number of chains to run}
  \item{nIter}{number of iterations per chain (not including burn-in and after thinning)}
  \item{nBurnin}{length of burn in, i.e. number of iterations to discard at the beginning.}
  \item{nThin}{thinning rate.  Must be a positive integer.
    Set \code{nThin} > 1 to save memory and computation time if \code{nIter} is large.
    Default is \code{max(1, floor(nIter / 10000))}
    which will only thin if there are at least 20000 simulations.}
  \item{DIC}{logical; if \code{TRUE} (default), obtain deviance samples and DIC estimates.
    The calculation of DIC done by WinBUGS; see the documentation of that package
    for more information.}
  \item{workingDirectory}{sets working directory during execution of this function;
    WinBUGS' or OpenBUGS' in- and output will be stored in this directory;
    if \code{NULL}, the current working directory is chosen.}
  \item{bugsDirectory}{directory that contains the WinBUGS executable, if engine = "WinBUGS";
    otherwise, this argument is ignored.}
  \item{debug}{if \code{FALSE} (default) and engine = "WinBUGS", WinBUGS is closed automatically
  	when the script has finished running, otherwise WinBUGS remains open for further investigation.}
  \item{engine}{which software engine to use to generate the posterior samples.  Options are
  	"WinBUGS" and "OpenBUGS".}
  \item{overRelax}{If \code{engine = OpenBUGS}, this allows "over-relaxation" of the sampling.
    See the WinBUGS documentation for more information, since this option is available through the 
    WinBUGS GUI in the model update menu.
    This option affects some but not all sampling methods in OpenBUGS.  }
  
}
\description{
	Obtains samples from the posterior distribution of the set of parameters, given the model and the data set.  
	The model specifies both the prior distributions of the parameters and the likelihood.  
}
\details{
  The samples
	are drawn using Markov Chain Monte Carlo, implemented in the packages WinBUGS and OpenBUGS.  The random seed in 
	BUGS is always the same, so with the same inputs the outputs will be exactly the same.  In order
	to perturb the output results, the initial values can be perturbed.
}

\value{
	An object of class "posterior" that contains the posterior samples.  
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
beta.prec <- 1e-6
# the following hyperparameter values result in a 
# very diffuse prior for the residual variance, one
# that has infinite mean and variance.
sigma2df <- 1  
sigma2scale <- 1


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
stackPost <- posteriorSamples (
  data = stackModelData, model = bugsModel, 
  inits = initialValues, nChains = nChains,
  parametersToSave = parameters.to.save,
  nIter = 10000, nBurnin = 1000, nThin = 5,
  engine = "WinBUGS")

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
# density estimates.
summary(stackPost)
densplot(stackPost)

# Obtain a "Bayesian p-value" for the coefficient for
# Air.Flow (here the Bayesian p-value means the 
# posterior probability that this coefficient is 
# negative).  The estimated Bayesian p-value is very 
# small or zero (varies by MCMC run), so the data 
# indicates that this coefficient is greater than 0.
betaSamples <- 
  getSamples(stackPost,"beta[2]")
sum(betaSamples < 0) / length(betaSamples)

}
