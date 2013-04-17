\name{seeds}
\alias{seeds}
\docType{data}
\title{Seed germination GLMM analysis}
\description{Seed germination data, as given in Crowder (1978). 
Here we fit a Bayesian overdispersed binomial model, assess convergence, and 
summarize the posterior distributions of the parameters.
  }
\usage{data(seeds)}
\format{
  A list with the following components
  \describe{
    \item{n}{A vector of the total number of seeds on each plate}
    \item{r}{A vector of the number of seeds that germinated on each plate}
    \item{x1}{A vector of the seed type for each plate}
    \item{x2}{A vector of the type of root extract for each plate}
  }
}
\source{
M. J. Crowder (1978). Beta-binomial ANOVA for proportions.
Applied Statistics, 27, 34-37.

N. E. Breslow and D. G. Clayton (1993).  Approximate Inference in Generalized Linear Mixed Models.  J. of 
the American Statistical Association, 88, 9-25.
}
\keyword{datasets}

\examples{

# Load the data
data(seeds)

# Specify a prior that is fairly diffuse relative 
# to this context
seedsPrior <- bhbm.prior ( 
  sigma2 = bayes.invChisq(df = 1, sigma0.sq = 1),
  fixed.coef = bayes.normal(mean = rep(0, 4), 
  cov = diag( rep(100, 4) ) ),
  common.glm = 2 )

# Specify the sampling parameters  
seedsSampler <- bhbm.sampler( nSamples = 1000,
  nThin = 10 )

# fit the model
seedsPost <- bhbm( data = seeds, 
  trials.formula = ~n, fixed.formula = r~x1+x2+x1*x2,
  overdispersion = "logit-normal", prior = seedsPrior,
  sampler = seedsSampler )

# perform convergence diagnostics
autocorr.plot(seedsPost)
traceplot(seedsPost)
geweke.diag(seedsPost)

# The convergence diagnostics look good, so 
# summarize the inferences
summary(seedsPost)

# The posterior coefficient estimates here are close
# to the maximum likelihood estimates given in 
# Breslow and Clayton (1993) in Table 3.
}
