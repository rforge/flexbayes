\name{bayes.uniformShrinkage}
\alias{bayes.uniformShrinkage}
\title{
  Uniform Shrinkage Prior Distribution
}
\description{
Set a non-informative uniform shrinkage prior distribution on a variance/scale parameter of a Bayesian model
 (see references below).
}
\usage{
bayes.uniformShrinkage(median.value=1.)
}
\arguments{
\item{median.value}{
a positive scalar specifying the median of the distribution.
}
}
\value{
an object of class \code{bayes.distribution} with the following slots:

\item{name}{a string equal to 
\code{"uniform shrinkage"}.}


\item{parameters}{a list containing:

 \bold{median:} the median of the distribution.
}
}
\details{
The uniform shrinkage distribution is one of the
 non-informative priors used for the variance parameters of a linear hierarchical model. Its density is given by
\code{f( tau^2 ) = tau0^2/( tau0^2 + tau^2 )^2 )}, where
\code{ tau0^2 } is the
\code{median.value} parameter.
}
\section{References}{
Christiansen, C. L. and Morris, C. N. (1997), "Hierarchical Poisson regression modeling", \emph{J. of the American Statistical Association}, 92: 618-632.

Daniels, M. J. (1999), "A prior for the variance in hierarchical models", \emph{The Canadian Journal of Statistics}, 27: 567-578.
}
\seealso{
\code{\link{bhlm.prior}},
\code{\link{bayes.invChisq}},
\code{\link{bayes.invWishart}},
\code{\link{bayes.massPoint}},
\code{\link{bayes.nonInfoPower}},
\code{\link{bayes.duMouchel}}.
}
\examples{
my.prior <- bayes.uniformShrinkage( 1.0 )
}
\keyword{Bayes}
% docclass is function
% Converted by Sd2Rd version 37351.
