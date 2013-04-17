\name{bayes.nonInfoPower}
\alias{bayes.nonInfoPower}
\title{
  Non-Informative Power Law Prior
}
\description{
Set a non-informative power law prior distribution on a variance/squared scale parameter of a Bayesian model.
In this case the distribution of the parameter is proportional to a (possibly zero)
 power of it.
}
\usage{
bayes.nonInfoPower(power=-1.)
}
\arguments{
\item{power}{
a scalar specifying the power to which the parameter distribution is proportional.
}
}
\value{
an object of class \code{bayes.distribution} with the following slots:

\item{name}{a string equal to 
\code{"non-informative power"}.}


\item{parameters}{a list containing: 

\bold{power}: the scalar specifying the power to which the
 parameter distribution is proportional.}

}
\section{References}{
Hobert, J. P., and G. Casella (1996). The effect of improper priors on Gibbs sampling in hierarchical linear mixed models. \emph{J. American Statistical Association}, 91: 1461-1473.
}
\seealso{
\code{\link{bhlm.prior}},
\code{\link{bayes.invChisq}},
\code{\link{bayes.invWishart}},
\code{\link{bayes.massPoint}},
\code{\link{bayes.duMouchel}},
\code{\link{bayes.uniformShrinkage}}.
}
\examples{
#set the prior distribution of sigma proportional to (sigma)^(-1)
sigma.prior <- bayes.nonInfoPower( -1.0 )

#set the prior to the Haar measure
sigma.prior <- bayes.nonInfoPower( -0.5 )
}
\keyword{Bayes}
% docclass is function
% Converted by Sd2Rd version 37351.
