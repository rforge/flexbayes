\name{nchain.posterior}
\alias{nchain}
\alias{nchain.posterior}
\alias{niter}
\alias{niter.posterior}
\alias{nvar}
\alias{nvar.posterior}
\title{Dimensions of posterior objects}

\usage{niter(x)
nvar(x)
nchain(x)}

\arguments{
\item{x}{A \code{posterior} object}
}

\description{
	\code{niter} returns the number of iterations in each chain.  \code{nvar} returns the 
	number of parameters (variables).  \code{nchain} returns the number of chains.
}
\seealso{
   \code{\link{start.posterior}},
   \code{\link{thin.posterior}},
   \code{\link{subset.posterior}},
   \code{\link{varnames.posterior}}.
}
