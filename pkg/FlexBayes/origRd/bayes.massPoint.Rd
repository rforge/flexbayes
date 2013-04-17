\name{bayes.massPoint}
\alias{bayes.massPoint}
\title{
  Mass Point Prior Distribution
}
\description{
Set an degenerate mass point prior distribution on a variance/scale parameter of a Bayesian model.
The distribution is concentrated in a single point, and hence makes the parameter a constant.
}
\usage{
bayes.massPoint(value=1)
}
\arguments{
\item{value}{
a scalar specifying the point where the distribution is concentrated.

 When specifying several mass point priors at the same time, such as in
the case of a hierarchical linear model with several groups,
the mass points associated to the different distributions of the parameters 
can be specified in a vector of length the number of parameters.
}
}
\value{
an object of class \code{bayes.distribution} with the following slots:

\item{name}{a string equal to \code{"mass point"}.}

\item{parameters}{a list containing:

 \bold{value:} the scalar specifying the point where the distribution is concentrated.
 }
}
\seealso{
\code{\link{bhlm.prior}},
\code{\link{bayes.invChisq}},
\code{\link{bayes.invWishart}},
\code{\link{bayes.nonInfoPower}},
\code{\link{bayes.duMouchel}},
\code{\link{bayes.uniformShrinkage}}.
}
\examples{
#give a constant value of 3.0 to the parameter.
my.massPoint <- bayes.massPoint( 3.0 )
}
\keyword{Bayes}
% docclass is function
% Converted by Sd2Rd version 37351.
