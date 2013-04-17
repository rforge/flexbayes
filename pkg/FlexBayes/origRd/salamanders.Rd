\name{salamanders}
\alias{salamanders}
\docType{data}
\title{Salamander mating data}
\description{Data on the mating success between two populations of salamanders, the roughbutts and 
the whitesides, from McCullogh and Nelder (1989).  This data is also analyzed by a number of other
authors, including Breslow and Clayton (1993).
  }
\usage{data(salamanders)}
\format{
  A list with the following components
  \describe{
    \item{x}{Indicators for the mating type: R/R, R/W, W/R, or W/W (the population for both the 
    female and male salamander)}
    \item{y}{The mating outcome for each of 3 experiments}
  }
}
\source{
P. McCullagh and J. A. Nelder (1989).  Generalized Linear Models, 2nd edition.  London: Chapman and Hall.

N. E. Breslow and D. G. Clayton (1993).  Approximate Inference in Generalized Linear Mixed Models.  J. of 
the American Statistical Association, 88, 9-25.
}
\keyword{datasets}

\examples{


# Specify the prior
salPrior <- bhbm.prior ( 
  sigma2 = bayes.invChisq(df = 2, sigma0.sq = 10),
  fixed.coef = bayes.normal( mean = rep(0, 4), 
  cov = rep(10, 4) ),
  random.var = bayes.invChisq(df = 2, sigma0.sq=10),
  common.glm = 2 )

# Specify the sampling parameters  
salSampler <- bhbm.sampler( nSamples = 1000,
  nThin = 10 )

# Load the data
data(salamanders)
# Prepare the data 
# Just analyze the first experiment
nMates <- dim(salamanders$y)[1]
salBhbm <- data.frame( y = salamanders$y[,1],
  x = salamanders$x, z = salamanders$z,
  trials = rep( 1, nMates ) )

# fit the model
fixedEff <- 
  y ~ x.R.R + x.R.W + x.W.R + x.W.W - 1
randomEff <- ~ z.fR1 + z.fR2 + z.fR3 + z.fR4 + 
  z.fR5 + z.fR6 + z.fR7 + z.fR8 + z.fR9 + z.fR10 + 
  z.fW1 + z.fW2 + z.fW3 + z.fW4 + z.fW5 + z.fW6 + 
  z.fW7 + z.fW8 + z.fW9 + z.fW10 + z.mR1 + z.mR2 +
  z.mR3 + z.mR4 + z.mR5 + z.mR6 + z.mR7 + z.mR8 +
  z.mR9 + z.mR10 + z.mW1 + z.mW2 + z.mW3 + z.mW4 +
  z.mW5 + z.mW6 + z.mW7 + z.mW8 + z.mW9 + z.mW10 - 1
# THIS TAKES A FEW MINUTES TO RUN
salPost <- bhbm( data = salBhbm, 
  fixed.formula = fixedEff, 
  random.formula = randomEff,
  prior = salPrior,
  sampler = salSampler, 
  trials.formula = ~trials,
  overdispersion = "logit-normal" )

# perform convergence diagnostics.  Note that by 
# default only the first thirty parameters are 
# displayed.
autocorr.plot(salPost)
effectiveSize(salPost)

# summarize the inferences.  Again only the first
# thirty parameters are displayed.  
summary(salPost)

}
