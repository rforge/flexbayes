\name{bhlm.sampler}
\alias{bhlm.sampler}
\title{
  Control and Initial Value Parameters for MCMC
}
\description{
Specifies
the number of parallel chains to be simulated, the control parameters 
governing each chain,
and the initial values of the model parameters for each chain. The resulting list must be
 passed to a call of \code{bhlm()} in order
to fit a Bayesian hierarchical linear model.
}
\usage{
bhlm.sampler(nBurnin = 1000, nSamples = 1000, 
  nThin = 1, nChains=1, init.point=
   "prior", params.init=NULL)
}
\arguments{
\item{nSamples}{number of samples drawn per chain (not including burn-in and after thinning)}
\item{nBurnin}{length of burn in, i.e. number of iterations to discard at the beginning of each chain.}
\item{nThin}{thinning rate.  Must be a positive integer.
    Set \code{nThin} > 1 to reduce the autocorrelation of the samples.}
\item{nChains }{
an integer specifying the number of parallel (independent) MCMC simulations (chains) to run.
In practice several chains should be run in order to test whether or not the chains have
converged to their limiting distribution. This test can be run
by calling \code{gelman.diag} on the output of the \code{bhlm} function.
}
\item{init.point }{
an string specifying the type of method to generate initial values for the 
\code{nChains}
simulations.

 \bold{"prior":} initial values for the fixed effect and second-level coefficients are drawn from t-distributions with 
location and scale given by the
prior specification for the coefficients; initial values for the outcome variance and random effect variances (or covariance matrix) are sampled from an appropriate inverse chi-squared (or inverse Wishart) distribution,
depending on the prior specification.

 \bold{"user's choice":} initial values for the coefficients and the variance components
 are explicitly provided
in \code{params.init}.
}
\item{params.init }{
a list specifying initial values for the MCMC simulations.
\code{params.init} is required only if 
\code{init.point} = "user's choice".
If  \code{nChains} > 1, then
\code{params.init}  must be a list of lists
containing the initial values for the parameters in each chain.
Otherwise it can simply be a list containing the initial values for the parameters.
Each chain's list must include the following elements:

 \bold{error.var:}
a vector, matrix or list specifying the initial values for the error variance parameter 
\bold{sigma^2}.
The type of input depends on the value of \code{common.error.var} in the call to \code{bhlm.prior} and the value of \code{nChains}. 
\code{error.var} must be a list if \code{common.error.var} = 0.

 \bold{random.coef:}
a vector, matrix or list specifying the initial values for the random effects coefficients
 (\bold{beta}).
If \code{random.coef} is a matrix or list,
 then the number of columns (or length, if a list)
must coincide with the number of groups in the model. The dimension (or number of rows) must
be the dimension of \bold{beta}.  If there are
no random effects in the model, these initial values will be ignored.

 \bold{fixed.coef: }
a vector specifying the initial values for the fixed effects coefficients (\bold{gamma}).  If there are
no fixed effects in the model, these initial values will be ignored.

 \bold{level2.coef:}
a vector specifying the initial values for the second level effects coefficients (\bold{alpha}).  If there are
no second-level effects in the model, these initial values will be ignored.

 \bold{random.var:}
a vector or matrix
specifying the initial values for the variances \bold{tau^2_k} or covariance matrix \bold{V} of the
random effects.  If there are
no random effects in the model, these initial values will be ignored.
}
}
\value{
a list containing the control parameters and initial values for the Gibbs sampler. 
}
\seealso{
\code{\link{bhlm}},
\code{\link{bhlm.prior}},
\code{\link{bhlm.likelihood}}.
}
\examples{

ortho <- as.data.frame( 
  Orthodont[ Orthodont$Sex == "Female", ] )

#the likelihood
ortho.lkhd <- bhlm.likelihood( type = "normal" )

#the priors
error.var.nu <- 3
error.var.sigma02 <- 0.5 * 0.5
error.var <- bayes.invChisq( df = error.var.nu, 
  sigma0.sq = error.var.sigma02 ) 

alpha.mean <- zero
alpha.Cov <- identity
level2.coef <- bayes.normal( 
  mean.vector = alpha.mean, covmat = alpha.Cov )

random.var.nu <- 3
random.var.sigma02 <- 1
random.var <- bayes.invChisq( df = random.var.nu, 
  sigma0.sq = random.var.sigma02 ) 

ortho.prior <- bhlm.prior( 
  error.var = error.var, random.var = random.var,
  level2.coef = level2.coef )

#the sampler parameters
ortho.sampler <- bhlm.sampler( init.point = "prior" )

#the call to bhlm
ortho.bhlm <- bhlm( 
  random.formula = distance ~ I(age - 11), 
  level2.formula = ~ 1,
  group.formula = ~ Subject, data = ortho, 
  prior = ortho.prior, 
  likelihood = ortho.lkhd, sampler = ortho.sampler )

traceplot(ortho.bhlm)
summary(ortho.bhlm)
}
\keyword{Bayes}
% docclass is function
% Converted by Sd2Rd version 37351.
