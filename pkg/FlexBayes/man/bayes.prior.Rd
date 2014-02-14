\name{bayes.prior}
\alias{bayes.normal}
\alias{bayes.t}
\alias{bayes.nonInformative}
\alias{bayes.beta}
\alias{bayes.gamma}
\alias{bayes.uniformShrinkage}
\alias{bayes.duMouchel}
\alias{bayes.nonInfoPower}
\alias{bayes.massPoint}
\alias{bayes.wishart}
\alias{bayes.invWishart}
\alias{bayes.invChisq}

\title{FlexBayes Prior Distribution Specification}

\description{
Specify the prior distribution for a FlexBayes model.
}

\usage{
  bayes.normal(mu = zero, sigma = identity, k0 = 1.0, ...)
  bayes.t(mu = zero, sigma = identity, df = 3, ...)
  bayes.nonInformative(mu = 0, sigma = 1, ...)
  bayes.beta(shape1, shape2, ...)
  bayes.gamma(shape, rate = 1.0, scale = 1.0/rate, ...)
  bayes.uniformShrinkage(median = 1.0)
  bayes.duMouchel(dispersion = 1.0)
  bayes.nonInfoPower(power = -1.0)
  bayes.massPoint(value = 1.0)
  bayes.wishart(df = 1, scale = 1.0)
  bayes.invWishart(df = 1, scale = 1.0)
  bayes.invChisq(df = 3, sigma0.sq = 1)
}

\arguments{
\item{mu}{a numeric vector of length \eqn{p} specifying the prior mean of the distribution. The default, \code{mu = zero}, sets the prior mean to be the \emph{zero} vector.}

\item{sigma}{a symmetric positive definite numeric matrix specifying the prior covariance structure of the distribution. Alternatively, if \code{sigma} is a single positive numeric value \eqn{\sigma^{2}} then the prior covariance matrix is taken to be \eqn{\sigma^{2} I_{p}} where \eqn{I_{p}} is the \eqn{p \times p} identity matrix. If \code{sigma} is a numeric vector containing \eqn{p} strictly positive values then the prior covariance matrix is taken to be a diagonal matrix with \code{sigma} along the main diagonal.  The default, \code{sigma = identity}, sets the prior covariance matrix to \eqn{I_{p}}.}

\item{k0}{a positive integer value specifying the \emph{equivalent} number of observations associated to the prior.  The larger \code{k0}, the larger the influence of the prior in the posterior distribution of the parameter. The effect of \code{k0} is to reduce the prior covariance matrix to \code{sigma / k0}.}

\item{df}{a positive integer value specifying the degrees of freedom.}

\item{shape1}{a positive numeric value specifying the first shape parameter for the beta distribution.}

\item{shape2}{a positive numeric value specifying the second shape parameter for the beta distribution.}

\item{shape}{a positive numeric value specifying the shape parameter for the gamma distribution.}

\item{rate}{a positive numeic value specifying the rate parameter.}

\item{scale}{a positive numeic value specifying the scale parameter.}

\item{median}{a positive numeric value.}

\item{dispersion}{a positive numeric value.}

\item{power}{a single numeric value.}

\item{value}{a single numeric value.}

\item{sigma0.sq}{a positive numeric value.}
}

\value{
a list with class \code{FBprior} containing the following elements:

\item{name}{the character string \code{"norm"}.}

\item{name}{the character string \code{"normal"}.}

\item{parameters}{a list containing the parameters of the prior distribution.}
}

\examples{
beta.prior <- bayes.normal(c(0.5, -0.5), covmat = identity)
}

\keyword{FlexBayes}
\keyword{prior}


