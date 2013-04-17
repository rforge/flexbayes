\name{geweke.diag}
\alias{geweke.diag}
\title{Geweke's convergence diagnostic}

\usage{
geweke.diag(x, maxVars = 30, \dots)
} 

\arguments{
  \item{x}{A \code{posterior} object}
	\item{maxVars}{Maximum number of parameters to display}
  \item{...}{Additional arguments for \code{geweke.diag} in the \code{coda} package}
}

\description{
Returns Z-scores for Geweke's diagnostic. Computation is done by calling 
the \code{geweke.diag} function in the \code{coda} library;
see the help file for \code{geweke.diag} in the \code{coda} package
for details and more optional arguments.  This can be found in the Available Help section of the Help menu when 
	\code{coda} is loaded, which occurs when \code{geweke.diag} is called 
	on an object of class \code{posterior}.
}

\references{
Geweke, J. (1993).  Evaluating the accuracy of sampling-based approaches
to calculation of posterior moments. In \emph{Bayesian Statistics 4}, eds.
 J. M. Bernardo, J. O. Berger, A. P. Dawid and A. F. M. Smith. 
Oxford: Oxford University Press.
}

\seealso{
   \code{\link{geweke.plot}},
   \code{\link{gelman.diag}},
   \code{\link{gelman.plot}}.
}
