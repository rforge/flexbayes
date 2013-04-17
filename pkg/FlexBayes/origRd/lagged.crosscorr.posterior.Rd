\name{lagged.crosscorr.posterior}
\alias{lagged.crosscorr}
\alias{lagged.crosscorr.posterior}
\title{Calculates the lagged cross-correlations for a posterior object}

\usage{lagged.crosscorr(x, maxVars = 30, lags = c(0, 1, 5, 10, 50), ...)}

\arguments{
\item{x}{A \code{posterior} object}
\item{maxVars}{Maximum number of parameters to display}
\item{lags}{The lag values at which to calculate the lagged cross-correlations}
\item{...}{Additional arguments to the \code{autocorr} function in the \code{coda} package}
}

\description{
	Calculates the lagged correlations of all pairs of parameters in a \code{posterior} object.  
	If there are multiple chains then the samples from all the chains are pooled.  The calculation is done by 
	calling \code{autocorr} in the \code{coda} package; see the \code{autocorr} section of the 
	\code{coda} help for details and 
	optional arguments. This can be found in the Available Help section of the Help menu when 
	\code{coda} is loaded, which occurs when \code{lagged.crosscorr} is called 
	on an object of class \code{posterior}.
}
\seealso{
   \code{\link{autocorr}},
   \code{\link{crosscorr}},
   \code{\link{crosscorr.plot}}.
}

