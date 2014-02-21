\name{fbprior}
\alias{fbprior}

\title{FlexBayes Prior Distribution}

\description{
Specify prior distributions for FlexBayes models.
}

\usage{
  fbprior(dstn, ...)
}

\arguments{
  \item{dstn}{a character string naimg the distribution.}

  \item{\dots}{parameters for the specified distribution.}
}

\value{
a list containing the following elements:

\item{name}{a character string namimg the distribution.}

\item{parameters}{a list containing the parameters of the distribution.}
}

\examples{
beta.prior <- fbprior("normal", mean = rep(0, 4), diag(4))
}

\keyword{FlexBayes}
\keyword{prior}


