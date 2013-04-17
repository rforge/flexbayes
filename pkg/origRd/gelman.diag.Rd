\name{gelman.diag}
\alias{gelman.diag}
\title{Gelman and Rubin's convergence diagnostic}

\usage{
gelman.diag(x, maxVars = 30, \dots)
} 

\arguments{
  \item{x}{A \code{posterior} object}
	\item{maxVars}{Maximum number of parameters for which to calculate the
	  diagnostic}
  \item{...}{Additional arguments for \code{gelman.diag} in the \code{coda} package}
}

\description{
Returns a 'potential scale reduction factor' along with associated confidence 
intervals.  Computation is done by calling the \code{gelman.diag} function in the \code{coda} library.
See the references given below or the help file for \code{gelman.diag} in the \code{coda} package
for details and more optional arguments.  The latter can be found in the Available Help section of the Help menu when 
	\code{coda} is loaded, which occurs when \code{gelman.diag} is called 
	on an object of class \code{posterior}.
}

\references{
Gelman, A. and D. B. Rubin (1992). Inference from iterative simulation
using multiple sequences. \emph{Statistical Science}, 7, 457-472.

}

\seealso{
   \code{\link{gelman.plot}},
   \code{\link{geweke.diag}},
   \code{\link{geweke.plot}}.
}
