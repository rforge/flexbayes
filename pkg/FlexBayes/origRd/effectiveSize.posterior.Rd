\name{effectiveSize.posterior}
\alias{effectiveSize}
\alias{effectiveSize.posterior}
\title{Calculate the effective sample sizes of the parameter posterior samples}

\usage{effectiveSize(x, maxVars = 30)}

\arguments{
\item{x}{A \code{posterior} object}
\item{maxVars}{Maximum number of parameters to display}
}

\description{
Calculates the effective sample size of the posterior samples for each parameter, which is the number of samples from
the chain
adjusted for autocorrelation of the chain.  If \code{x} has multiple chains the effective sample sizes are summed across chains.  
The calculation is done by calling the \code{effectiveSize} function in the \code{coda} package.  For details on the calculation, 
	see the \code{effectiveSize} section of the \code{coda} 
	help.  This is available in the Available Help section of the Help menu when 
	\code{coda} is loaded, which occurs when \code{effectiveSize} is called on an object of class \code{posterior}.
}
\seealso{
\code{\link{autocorr}},
\code{\link{autocorr.plot}}.
}

