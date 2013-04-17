\name{autocorr.posterior}
\alias{autocorr}
\alias{autocorr.posterior}
\title{Calculate autocorrelations of a posterior object}

\usage{autocorr(x, maxVars = 30, lags = c(0, 1, 5, 10, 50), ...)}

\arguments{
\item{x}{A \code{posterior} object}
\item{maxVars}{Maximum number of parameters for which the autocorrelation is calculated}
\item{lags}{The lag lengths at which autocorrelations should be calculated}
\item{...}{Other arguments to \code{autocorr.diag} in the \code{coda} package}
}

\description{
	For each parameter in \code{x}, calculates the lagged autocorrelations of the posterior
	samples.  If \code{x} has 
	multiple chains, then the autocorrelations are averaged over the chains.
	The computation is done by calling \code{autocorr.diag} in the \code{coda} package. For 
	more information on the 
	optional arguments, see the \code{autocorr.diag} section
	of the \code{coda} 
	help.  This is available in the Available Help section of the Help menu when 
	\code{coda} is loaded, which occurs when \code{autocorr} is called 
	on an object of class \code{posterior}.
}
\seealso{
\code{\link{crosscorr}},
\code{\link{autocorr.plot}}.
}

