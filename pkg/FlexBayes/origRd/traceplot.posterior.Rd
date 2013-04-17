\name{traceplot.posterior}
\alias{traceplot}
\alias{traceplot.posterior}
\title{Trace plot of the posterior samples}

\usage{traceplot(x, maxVars = 30, ...)}

\arguments{
  \item{x}{A \code{posterior} object}
	\item{maxVars}{Maximum number of parameters to display}
  \item{...}{Additional arguments for \code{traceplot} in the \code{coda} package}
}

\description{
Creates a plot of the sampled value versus the iteration for each parameter. These plots 
are called "trace plots", or alternatively "time series plots". 
When there are multiple chains, their traces are superimposed.  The plots are generated
 by calling \code{traceplot} in the \code{coda} package; see 
the \code{traceplot} section of the \code{coda} help
  for more optional arguments.  This can be found in the Available Help section of the Help menu when 
	\code{coda} is loaded, which occurs when \code{traceplot} is called 
	on an object of class \code{posterior}.
}
\seealso{
   \code{\link{autocorr.plot}},
   \code{\link{geweke.plot}},
   \code{\link{gelman.plot}}.
}

