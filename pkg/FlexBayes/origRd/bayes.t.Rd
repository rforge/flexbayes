\name{bayes.t}
\alias{bayes.t}
\title{
  Student-t Prior Distribution
}
\description{
Specify a multivariate t prior distribution for a parameter of a Bayesian model.
}
\usage{
bayes.t(mean.vector = zero, covmat = identity, 
  df = 3)
}
\arguments{
\item{mean.vector }{
a vector specifying the location of the distribution.
If the prior location is the zero vector, one can literally
set \code{mean.vector = zero}.
}
\item{covmat}{
the scale structure of the t distribution.

 If \code{covmat} is specified as a positive scalar, then
the scale matrix is assumed to be of the form the scalar times the identity
matrix.

 If \code{covmat} is specified as a positive vector, then
the scale matrix is assumed to be diagonal with main diagonal given by
the vector \code{covmat}.

If \code{covmat} is a symmetric positive-definite matrix then this 
is taken to be the scale matrix of the t distribution.

If the prior scale matrix is the identity, one can literally set
\code{covmat = identity}.
}
\item{df}{
the prior degrees of freedom.
}
}
\value{
an object of class \code{bayes.distribution} with the following slots:

\item{name}{a string equal to \code{"t"}}. 

\item{parameters}{a list containing:

 \bold{mean.vector}: the location of the t distribution.

 \bold{covariance matrix:} the scale matrix of the t distribution.

 \bold{t degrees of freedom:} the degrees of freedom.}

}
\section{References}{
Gelman, A., J. B. Carlin, H. S. Stern, and D. B. Rubin (2004). \emph{ Bayesian Data Analysis}, 2nd Edition. Boca Raton, Florida: Chapman & Hall.
}
\seealso{
\code{\link{bhlm.prior}}, \code{\link{bayes.normal}}.
}
\examples{
my.t <- bayes.t( df = 5 )

#two parameter vectors, different priot locations, same prior scale.
location <- matrix( c( -1, 0, 2, 1, 1, 2 ), ncol = 2 )
my.t <- bayes.t( mean.vector = location, covmat = 0.25 * diag( 3 ), df = 3 )
}
\keyword{Bayes}
% docclass is function
% Converted by Sd2Rd version 37351.
