\name{bayes.normal}
\alias{bayes.normal}

\title{Multivariate Normal Prior Distribution Specification}

\description{
Specify a multivariate normal distribution as a prior in a FlexBayes model.
}

\usage{
bayes.normal(mu = zero, sigma = identity, k0 = 1.0, ...)
}

\arguments{
\item{mu}{a numeric vector of length \eqn{p} specifying the prior mean of the distribution. The default, \code{mu = zero}, sets the prior mean to be the \emph{zero} vector.}

\item{sigma}{a symmetric positive definite numeric matrix specifying the prior covariance structure of the distribution. Alternatively, if \code{sigma} is a single positive numeric value \eqn{\sigma^{2}} then the prior covariance matrix is taken to be \eqn{\sigma^{2} I_{p}} where \eqn{I_{p}} is the \eqn{p \times p} identity matrix. If \code{sigma} is a numeric vector containing \eqn{p} strictly positive values then the prior covariance matrix is taken to be a diagonal matrix with \code{sigma} along the main diagonal.  The default, \code{sigma = identity}, sets the prior covariance matrix to \eqn{I_{p}}.}

\item{k0}{a positive integer value specifying the \emph{equivalent} number of observations associated to the prior.  The larger \code{k0}, the larger the influence of the prior in
the posterior distribution of the parameter. The effect of \code{k0} is to reduce the prior covariance matrix to \code{sigma / k0}.}
}

\value{
a list with class \code{fbdstn} containing the following elements:

\item{name}{the character string \code{"norm"}.}

\item{name}{the character string \code{"normal"}.}

\item{parameters}{a list with elements:

\describe{
\item{mu:}{a numeric vector of length \eqn{p} containing the prior mean.}

\item{sigma:}{a numeric matrix containing the prior covariance matrix.}

\item{k0:}{a positive integer value giving the equivalent number of observations associated with the prior.}
}
}
}

\examples{
beta.prior <- bayes.normal(c(0.5, -0.5), covmat = identity)
}

\keyword{FlexBayes}
\keyword{prior}


