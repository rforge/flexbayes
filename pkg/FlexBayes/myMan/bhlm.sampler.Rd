\name{bhlm.sampler}
\alias{bhlm.sampler}

\title{MCMC Control and Initial Value Parameters}

\description{
  Specify the number of parallel chains to be simulated, the control parameters governing each chain, and the initial values of the model parameters for each chain.
}

\usage{
bhlm.sampler(nBurnin = 1000, nSamples = 1000, nThin = 1, nChains = 1,
             init.point = "prior", params.init = NULL)
}

\arguments{
  \item{nBurnin}{a positive integer value specifying the length of the burn in period (i.e., the number of iterations to discard from the beginning of each chain).}

  \item{nSamples}{a positive integer value specifying the number of samples drawn per chain (not including burn-in and after thinning).}

  \item{nThin}{a positive integer specifying the thinning rate.  Setting \code{nThin} to a value greater than 1 is highly likely to reduce the autocorrelation of the samples.}

  \item{nChains}{a positive integer value specifying the number of parallel (independent) MCMC simulations (chains) to run.  In practice several chains should be run in order to test whether the chains have converged to their limiting distribution.  This test can be run by calling \code{gelman.diag} on the output of the \code{bhlm} function.}

  \item{init.point}{a string specifying the method used to generate initial values for the simulations.  \describe{ \item{\code{"prior"}}{Initial values for the fixed effect and second-level coefficients are drawn from t distributions with location and scale given by the prior specification for the coefficients.  Initial values for the outcome variance and random effect variances (or covariance matrix) are sampled from an appropriate inverse chi-squared (or inverse Wishart) distribution, depending on the prior specification.} \item{\code{"user"}}{Initial values for the coefficients and the variance components are explicitly provided in the \code{params.init} argument to this function.}}}

  \item{params.init}{a list specifying initial values for the MCMC simulations to be used when \code{init.point = "user"}.  If \code{nChains} is greater than 1, then \code{params.init} must be a list of lists containing the initial values for the parameters in each chain.  Otherwise it can simply be a list containing the initial values for the parameters.  Each chain's list must include the following elements: \describe{\item{\code{error.var}}{a vector, matrix or list specifying the initial values for the error variance parameter \eqn{\sigma^{2}}.  The type of input depends on the value of \code{common.error.var} in the call to \code{bhlm.prior} and the value of \code{nChains}.  \code{error.var} must be a list if \code{common.error.var} = 0.}  \item{\code{random.coef}}{a vector, matrix or list specifying the initial values for the random effects coefficients \eqn{\beta}. If \code{random.coef} is a matrix or list, then the number of columns (or length, if a list) must coincide with the number of groups in the model. The dimension (or number of rows) must be the dimension of \eqn{\beta}.  If there are no random effects in the model, these initial values will be ignored.} \item{\code{fixed.coef}}{a vector specifying the initial values for the fixed effects coefficients \eqn{\gamma}.  If there are no fixed effects in the model, these initial values will be ignored.}  \item{\code{level2.coef}}{a vector specifying the initial values for the second level effects coefficients \eqn{\alpha}.  If there are no second-level effects in the model, these initial values will be ignored.}  \item{\code{random.var}}{a vector or matrix specifying the initial values for the variances \eqn{\tau_{k}^{2}} or covariance matrix \eqn{V} of the random effects.  If there are no random effects in the model, these initial values will be ignored.}}}
}

\value{
a list containing the control parameters and initial values for the Gibbs sampler. 
}

\seealso{
\code{\link{bhlm}}
}

\keyword{FlexBayes}
\keyword{prior}
