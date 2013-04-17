\name{preterm}
\alias{preterm}
\docType{data}
\title{Relationship of preterm birth to maternal DDE levels}
\description{Data on the relationship between maternal blood levels of DDE (a DDT metabolite)
and preterm birth, from Longnecker et al. (2001).  In the example below we fit a logistic regression model 
to the data to assess the relationship of DDE to preterm birth, while accounting for potential confounding
factors such as maternal smoking.
  }
\usage{data(preterm)}
\format{
  A list with the following components
  \describe{
    \item{x}{The maternal DDE level}
    \item{y}{An indicator for preterm birth}
    \item{z1}{Potential confounding factor (maternal smoking, etc.)}
    \item{z2}{Potential confounding factor}
    \item{z3}{Potential confounding factor}
    \item{z4}{Potential confounding factor}
    \item{z5}{Potential confounding factor}
  }
}
\source{
M. P. Longnecker, M. A. Klebanoff, H. Zhou and J. W. Brock (2001). Association between maternal serum
concentration of the DDT metabolite DDE and preterm and small-for-gestational-age babies at birth.
The Lancet, 358, 110-115.

}
\keyword{datasets}

\examples{

# Specify the prior
pretermPrior <- bhbm.prior ( 
  sigma2 = bayes.invChisq(df = 3, sigma0.sq = 1),
  fixed.coef = bayes.normal( mean = rep(0, 2), 
  cov = diag( rep(10, 2) ) ),
  common.glm = 2 )

# Specify the sampling control arguments, including 
# initial values for the parameters.  Notice that
# initial values are provided for the random 
# coefficients and random effect variance, although
# these parameters will not be in the model that
# is fit.  The initial values for these parameters
# will be ignored.
pretermSampler <- bhbm.sampler( nSamples = 1000,
  nThin = 1, nBurnin = 1000, 
  init.point = "user's choice",
  params.init = list( sigma2 = 2.9, 
  fixed.coef = c(-3.2, 0.02), random.coef = 0,
  random.var = 1 ) )

# Load the data
data(preterm)
pretermBhbm <- preterm
# Add the number of "trials" to the dataset, namely one for each birth
nBirths <- dim(preterm)[1]
pretermBhbm$trials <- rep(1, nBirths)

# fit the model
fixedEff <- 
  y ~ x
# THIS TAKES A FEW MINUTES TO RUN
pretermPost <- bhbm( data = pretermBhbm, 
  fixed.formula = fixedEff, 
  prior = pretermPrior, 
  trials.formula = ~ trials,
  sampler = pretermSampler, 
  overdispersion = "logit-normal" )

# perform convergence diagnostics.  Unfortunately
# there is a problem with the mixing here due
# to the large number of augmentation parameters
# (theta's), so we cannot obtain good posterior 
# estimates from this run.  A much longer run would
# be necessary; even 100x as many samples does not
# fix the problem.
traceplot(pretermPost)

}
