\name{bayes.invWishart}
\alias{bayes.invWishart}
\title{
  Inverse Wishart Prior Distribution
}
\description{
Set an inverse Wishart prior distribution on a covariance matrix in a Bayesian model.
The distribution is given by the scale matrix and degrees of freedom. It corresponds to the
generalization of the inverse chi-squared distribution to higher dimensions.
}
\usage{
bayes.invWishart(df=1, scale)
}
\arguments{
\item{scale}{
a symmetric positive-definite matrix specifying the scale matrix of the distribution.

 If \code{scale} is specified as a positive scalar, then
the scale matrix is assumed to be of the form the scalar times the identity
matrix.

 If \code{scale} is specified as a positive vector, then
the scale matrix is assumed to be diagonal with main diagonal given by
the vector \code{scale}.

If the prior scale is the identity matrix, one can literally set
\code{scale = identity}.
}
\item{df}{
the degrees of freedom associated to the distribution.
}
}
\value{
an object of class \code{bayes.distribution} with the following slots:

\item{name}{a string equal to \code{"invWishart"}.}

\item{parameters}{a list containing:

 \bold{degrees of freedom:} the degrees of freedom of the distribution.

 \bold{scale:} the scale matrix of the distribution.

}}
\section{References}{
Gelman, A., J. B. Carlin, H. S. Stern, and D. B. Rubin (2004). \emph{ Bayesian Data Analysis}, 2nd edition. 
Boca Raton, FL: Chapman & Hall.
}
\seealso{
\code{\link{bhlm.prior}},
\code{\link{bayes.invChisq}},
\code{\link{bayes.massPoint}},
\code{\link{bayes.nonInfoPower}},
\code{\link{bayes.duMouchel}},
\code{\link{bayes.uniformShrinkage}}.
}
\details{The distribution parameterization is as specified in Gelman, Carlin, Stern, and Rubin (2004).
}
\examples{
#a simple prior
scale.prior <- bayes.invWishart( df = 6, scale = diag( 6 )/6 )


#an empirical and informative prior
y.df <- 6
y.cov <- var( y )
y.inv <- solve( y.df * y.cov )
#symmetrized it
y.inv <- 0.5 * ( y.inv + t( y.inv ) )
#set the prior
y.cov.prior <- bayes.invWishart( df = y.df, scale = y.inv )
}
\keyword{Bayes}
% docclass is function
% Converted by Sd2Rd version 37351.
