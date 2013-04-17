\name{posterior}
\alias{posterior}
\title{Create an object of class posterior}

\usage{posterior(x)}

\arguments{
  \item{x}{An object of class mcmc, mcmc.list, matrix, or data.frame that
  	contains the posterior samples.}
}

\description{
Creates an object of class posterior.  \code{x} may be an mcmc or mcmc.list object (defined in
the coda package) or a matrix or data frame containing the samples.  If it is a matrix,
the columns should be the parameters, and the rows should be the iterations of the chain.
The columns should be named according to the parameter names.

If \code{x} is not of class mcmc.list, the resulting object will have a single chain.  
}

\value{
An object of class posterior, containing the posterior samples given in x.  
}

