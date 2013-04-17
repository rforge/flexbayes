\name{bhpm.sampler}
\alias{bhpm.sampler}
\title{
  Control and initial value parameters for Poisson model MCMC
}
\description{
Returns a list specifying
the number of parallel chains to be simulated, the control parameters 
governing each chain,
and the initial values of the model parameters for each chain. This list must be
 passed to a call of \code{\link{bhpm}} in order
to fit a Bayesian hierarchical Poisson model.
}
\usage{
bhpm.sampler(nBurnin = 1000, nSamples = 1000, 
  nThin = 1, nChains=1, update.cov=1, 
  init.point="prior", 
  params.init=NULL)
}
\arguments{
\item{nSamples}{number of samples drawn per chain (not including burn-in and after thinning)}
\item{nBurnin}{length of burn-in, i.e. the number of iterations to discard at the beginning.}
\item{nThin}{thinning rate.  Must be a positive integer.
    Set \code{nThin} > 1 to reduce the autocorrelation of the samples.}
\item{nChains }{
an integer specifying the number of parallel (independent) MCMC simulations (chains) to run.
It is common to run several chains in order to test whether or not the chains have
converged to their limit stationary distribution. This test can be run
by calling \code{gelman.diag} on the output of the \code{bhpm} function.
}
\item{update.cov }{
an integer specifying whether the proposal distribution covariance associated with the random and fixed effects
coefficients are to be updated at each iteration of the Metropolis-Hastings sampler. The default
is \code{update.cov = 1}. If \code{update.cov = 0}, then
the proposal covariances are fixed during the first iteration of the sampler, and are never
updated again during the simulation.
}
\item{init.point }{
a string specifying the type of method to generate initial values for the 
\code{nChains}
simulations.

 \bold{"prior":} initial values are drawn from a t-distribution with 
location and scale given by the
prior specification for the coefficients, from an inverse chi-squared (or inverse Wishart) 
based on the prior specification for the variances, and from a uniform shrinkage
distribution for the \bold{xi} or \bold{sigma^2} parameter 
(see \code{\link{bhpm.prior}} for details).

 \bold{"user's choice":} initial values for the coefficients and variances
 are explicitly provided
in \code{params.init}.
}
\item{params.init }{
a list specifying the initial values for the MCMC simulations.
\code{params.init} is required only if 
\code{init.point = "user's choice"}.
If  \code{nChains > 1}, then
\code{params.init}  must be a list of lists
containing the initial values for the parameters in each chain.
Each chain list must include all of the following elements (except where noted):

 \bold{xi:}
a scalar or vector specifying the initial values for the \bold{xi}
parameter (for the "gamma-conj" overdispersion case).
The type of input depends on the value of \code{common.glm} in the 
call to \code{\link{bhpm.prior}}. 

 \bold{sigma2:}
a scalar or vector specifying the initial values for the variance parameter \bold{sigma^2} 
(for the "log-normal" Poisson model).
The type of input depends on the value of \code{common.glm} in the call to \code{\link{bhpm.prior}}.  
Either \code{xi} or
\code{sigma2} (as appropriate for the 
model) must be specified in \code{params.init}.

 \bold{random.coef:}
a vector or matrix specifying the initial values for the random effects coefficients
 (\bold{beta}).
If \code{random.coef} is a matrix,
 then the number of columns
must coincide with the number of groups in the model, and the number of rows must
be the dimension of \bold{beta}. 
\code{random.coef} must be specified; if there are no random effects in the model then the value of \code{random.coef} will
be ignored.

 \bold{fixed.coef: }
a vector specifying the initial values for the fixed effects coefficients (\bold{gamma}).
Must be included; if there are no fixed effects in the model then the value of \code{fixed.coef} will
be ignored.

 \bold{level2.coef:}
a vector specifying the initial values for the second level effects coefficients (\bold{alpha}).  
Only needs to be included if there are second-level coefficients in the model to be fit.

 \bold{random.var:}
a scalar or vector (depending on the prior) 
specifying the initial value for the variance structure \bold{tau^2} of the
random effects coefficients.
}
}
\value{
a list containing the control parameters and initial values for the Metropolis-Hastings sampler.
}
\seealso{
\code{\link{bhpm}},
\code{\link{bhpm.prior}}.
}
\examples{
#fit the pump data
pump.params <- list( xi = 0.225, random.coef = 
  c( -1.0, -1.0 ), 
  fixed.coef = 0, random.var = c(1,1) )
pump.sampler <- bhpm.sampler( nBurnin = 1000, 
  nSamples = 1000, nThin = 50, update.cov = 1, 
  init.point = "user's choice", 
  params.init = pump.params )
  
xi.prior = bayes.uniformShrinkage (0.5)
random.var.prior = bayes.invWishart( df = 3, scale = diag( c(1,1) ) )

pump.prior = bhpm.prior ( xi = xi.prior,
  random.var =  random.var.prior,
  common.glm = 2 )

pump.exposure <- ~ e
pump.random <- z ~ 1 + x

##########################################
## call bhpm to fit a gamma-conjugate model
pump.bhpm <- bhpm( random.formula = pump.random, 
  exposure.formula = pump.exposure, data = pumps, 
  prior = pump.prior, sampler = pump.sampler,
  overdispersion = "gamma-conj" )

}
\keyword{Bayes}
