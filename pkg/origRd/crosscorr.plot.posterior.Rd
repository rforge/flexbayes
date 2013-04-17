\name{crosscorr.plot.posterior}
\alias{crosscorr.plot}
\alias{crosscorr.plot.posterior}
\title{Plot crosscorrelations of the posterior samples}

\usage{crosscorr.plot(x, maxVars = 30, ...)}

\arguments{
\item{x}{A \code{posterior} object}
\item{maxVars}{Maximum number of parameters to display}
\item{...}{Other arguments to \code{crosscorr.plot} in the \code{coda} package}
}

\description{
  Creates an image plot of the pairwise correlation matrix for the parameters in \code{x}.  If \code{x} has multiple
  chains, then the samples from all of the chains are pooled.
  Creates the plot by calling the \code{crosscorr.plot} function in the \code{coda}
  package.  For more information on the 
	optional arguments to this function, see the \code{crosscorr.plot} section
	of the \code{coda} 
	help.  This is available in the Available Help section of the Help menu when 
	\code{coda} is loaded, which occurs when \code{crosscorr.plot} is called on an object
	of class \code{posterior}.
}
\seealso{
\code{\link{crosscorr}},
\code{\link{autocorr.plot}}.
}

