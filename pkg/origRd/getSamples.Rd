\name{getSamples}
\alias{getSamples}
\title{Posterior Samples as a Matrix}

\usage{
  getSamples(post, param = NULL, chainIndex = NULL)
}
\arguments{
  \item{post}{An object of class \code{posterior}.}
  \item{param}{The names or indexes of the parameters for which to obtain the samples, 
  as a single value or a vector.
  If NULL, use all of the parameters.}
  \item{chainIndex}{The index of the chain for which to obtain the samples.  If NULL,
    obtain the samples from all of the chains concatenated to give a single matrix.}
}
\description{
	Returns the samples from a \code{posterior} object as a matrix with \code{nvar(post)} columns 
	and \code{niter(post) * nchain(post)} rows.
}

\value{
	A matrix that contains the posterior samples.  
}

