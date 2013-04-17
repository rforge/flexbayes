\name{geweke.plot}
\alias{geweke.plot}
\title{Geweke-Brooks plot}

\usage{
geweke.plot(x, maxVars = 6, \dots)
} 

\arguments{
  \item{x}{A \code{posterior} object}
	\item{maxVars}{Maximum number of parameters to display.  The default
	  is small because these diagnostics can be time-consuming to calculate.}
  \item{...}{Additional arguments for \code{geweke.plot} in the \code{coda} package}
}

\description{
Creates a plot of Geweke's Z-score.
Computation is done by calling the \code{geweke.plot} function in the \code{coda} library; 
see the help file for \code{geweke.plot} in the \code{coda} package
for details and more optional arguments.  This can be found in the Available Help section of the Help menu when 
	\code{coda} is loaded, which occurs when \code{geweke.plot} is called 
	on an object of class \code{posterior}.
}

\seealso{
   \code{\link{geweke.diag}},
   \code{\link{gelman.diag}},
   \code{\link{gelman.plot}}.
}
