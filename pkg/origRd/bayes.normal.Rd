\name{bayes.normal}
\alias{bayes.normal}
\title{
  Multivariate Normal Prior Distribution
}
\description{
Specify a multivariate normal prior distribution for a parameter of a Bayesian model.
}
\usage{
bayes.normal(mean.vector = zero, covmat= identity, 
  k0 = 1)
}
\arguments{
\item{mean.vector }{
a vector specifying the mean of the distribution.
If the prior mean is the zero vector, one can literally
set \code{mean.vector = zero}.
}
\item{covmat}{
the covariance structure of the
normal distribution.

 If \code{covmat} is specified as a positive scalar, then
the covariance matrix is assumed to be of the form the scalar times the identity
matrix.

 If \code{covmat} is specified as a positive vector, then
the covariance matrix is assumed to be diagonal with main diagonal given by
the vector \code{covmat}.

If \code{covmat} is a symmetric positive-definite matrix, then this is taken
to be the covariance matrix of the normal distribution.

If the prior covariance is the identity matrix, one can literally set
\code{covmat = identity}.
}
\item{k0}{
an integer specifying the \emph{equivalent} number of observations associated to the prior.
The larger \code{k0} the larger the influence of the prior in
the posterior distribution of the parameter. The effect of \code{k0} is
to reduce the prior covariance matrix to \code{covmat / k0}.
}
}
\value{
an object of class \code{bayes.distribution} with the following slots:

\item{name}{ a string equal to \code{"normal"}. }

\item{parameters}{a list containing:

\bold{mean.vector:} the mean of the normal distribution.

 \bold{covariance matrix:} the covariance of the normal distribution.

 \bold{k0:} the equivalent number of observations associated to the prior.}

}
\seealso{
\code{\link{bayes.t}}, \code{\link{bhlm.prior}}.
}
\section{References}{
Gelman, A., J. B. Carlin, H. S. Stern, and D. B. Rubin (2004). \emph{ Bayesian Data Analysis}, 2nd Edition. Boca Raton, Florida: Chapman & Hall.
}
\examples{
my.normal <- bayes.normal( mean.vector = c( 0.5, -0.5 ), covmat = identity )
}
\keyword{Bayes}
% docclass is function
% Converted by Sd2Rd version 37351.
