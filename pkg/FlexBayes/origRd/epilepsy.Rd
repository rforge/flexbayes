\name{epilepsy}
\alias{epilepsy}
\docType{data}
\title{Epilepsy drug hierarchical Poisson modeling}
\description{Data on the effect of a drug on the occurrence of seizures, as given in Thall and 
Vail (1990).  Here we fit a Bayesian hierarchical version of a Poisson-outcome model
fit to the data by Breslow and Clayton (1993) (their model III).  We assess convergence and 
summarize the posterior distributions of the parameters.
  }
\usage{data(epilepsy)}
\format{
  A list with the following components
  \describe{
    \item{y}{The seizure counts for each patient in the two weeks proceeding each of four visits}
    \item{visit4}{An indicator of visit four, for each visit}
    \item{baseline}{The number of seizures in the 8 weeks preceeding the trial, for each patient}
    \item{treatment}{An indicator of treatment (as opposed to control) for each patient}
    \item{age}{The age of each patient}
  }
}
\source{
P. F. Thall and S. C. Vail (1990). Some covariance models for longitudinal count data with overdispersion.
Biometrics, 46, 657-671.

N. E. Breslow and D. G. Clayton (1993).  Approximate Inference in Generalized Linear Mixed Models.  J. of 
the American Statistical Association, 88, 9-25.
}
\keyword{datasets}

\examples{


# Specify the prior
epilepsyPrior <- bhpm.prior ( 
  sigma2 = bayes.nonInfoPower(-1), 
  common.glm = 2)

# Specify the sampling parameters  
epilepsySampler <- bhpm.sampler( nSamples = 1000,
  nThin = 10, init.point = "user's choice",
  params.init = list( sigma2 = 1, 
  fixed.coef = rep(0, 1), random.var = 1,
  random.coef = 0, level2.coef = rep(0, 5) ) )

# Load the data
data(epilepsy)
# Prepare the data 
nSubj <- dim(epilepsy$y)[1]
base <- log ( epilepsy$baseline / 4 )
epilepsyBhpm <- data.frame( 
  y = as.vector( epilepsy$y ),
  subj = rep( (1:nSubj), each = 4 ),
  visit4 = rep( epilepsy$visit4, nSubj ),
  base = rep( base, each = 4 ),
  treat = rep( epilepsy$treatment, each = 4 ),
  logAge = rep( log( epilepsy$age ), each = 4 ) )  

# fit the model
fixedEff <- y ~ visit4 - 1
level2Eff <- ~ treat + base + logAge + treat * base
randomEff <- ~ 1
# THIS TAKES A FEW MINUTES TO RUN
epilepsyPost <- bhpm( data = epilepsyBhpm, 
  fixed.formula = fixedEff, 
  random.formula = randomEff,
  group.formula = ~ subj,
  level2.formula = level2Eff,
  prior = epilepsyPrior,
  sampler = epilepsySampler, 
  overdispersion = "log-normal" )


# perform convergence diagnostics.  Note that by 
# default only the first thirty parameters are 
# displayed.
autocorr.plot(epilepsyPost)
traceplot(epilepsyPost)
geweke.diag(epilepsyPost)

# The convergence diagnostics look good, so 
# summarize the inferences.  Again only the first
# thirty parameters are displayed.  The (Intercept):#
# parameters are the random effects.  
# The results are similar to those reported by 
# Breslow and Clayton (1993).
summary(epilepsyPost)
# Look at the names of the rest of the parameters in 
# the model.  The lambda parameters are data 
# augmentation parameters for the overdispersion
# (one lambda parameter for each data point).
varnames(epilepsyPost)



}
