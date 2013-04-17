\name{bayes.duMouchel}
\alias{bayes.duMouchel}
\title{
  Non-Informative DuMouchel prior
}
\description{
Set a non-informative prior distribution on a variance/scale parameter of a Bayesian model.
In this case the distribution of the parameter is given by a \emph{DuMouchel} prior
 (see references below).
}
\usage{
bayes.duMouchel(dispersion=1.)
}
\arguments{
\item{dispersion}{
a positive scalar characterizing the distribution.
}
}
\value{
an object of class \code{bayes.distribution}.
}
\details{
The DuMouchel distribution is one of the
 non-informative priors used on a linear hierarchical model. Its density is given by
\code{f( tau^2 ) = tau0/( 2 * tau * ( tau0 + tau )^2 )}, where
\code{ tau0 } is the
\code{dispersion} parameter.
}
\section{References}{
Daniels, M. J., "A prior for the variance in hierarchical models", \emph{The Canadian Journal of Statistics}, vol. 27, 3: 567-578, 1999.
}
\seealso{
\code{\link{bhlm.prior}},
\code{\link{bayes.invChisq}},
\code{\link{bayes.invWishart}},
\code{\link{bayes.massPoint}},
\code{\link{bayes.nonInfoPower}},
\code{\link{bayes.uniformShrinkage}}.
}
\examples{
my.prior <- bayes.duMouchel( 1.0 )
}
\keyword{Bayes}
