\name{bayes.invChisq}
\alias{bayes.invChisq}
\title{
  Inverse Chi-square Prior Distribution
}
\description{
Specify an inverse chi-square prior distribution for a variance / squared scale parameter of a 
Bayesian model. The distribution is given by the degrees of freedom and prior squared scale.
}
\usage{
bayes.invChisq(df=3, sigma0.sq=1)
}
\arguments{
\item{df}{
the degrees of freedom of the distribution.
}
\item{sigma0.sq }{
a positive scalar specifying the square of the scale of the distribution.
}
}
\value{
an object of class \code{bayes.distribution} with the following slots:

\item{name}{ a string equal to \code{"invChisq"}.}

\item{parameters}{ a list containing:

 \bold{df:} the degrees of freedom of the distribution.

 \bold{sigma0.sq:} the square of the scale of the distribution.}

}
\section{References}{
Gelman, A., J. B. Carlin, H. S. Stern, and D. B. Rubin (2004). \emph{ Bayesian Data Analysis}, 2nd Edition. 
Boca Raton, Florida: Chapman & Hall.
}
\seealso{
\code{\link{bhlm.prior}}
}
\details{This distribution is the scaled inverse chi-squared distribution as specified in Gelman,
Carlin, Stern, and Rubin (2004).}
\examples{
#setting a non-informative prior
my.scale <- bayes.invChisq()

#setting an informative prior
my.scale <- bayes.invChisq( df = 3, sigma0.sq = (0.025)^2 )
}
\keyword{Bayes}
% docclass is function
% Converted by Sd2Rd version 37351.
