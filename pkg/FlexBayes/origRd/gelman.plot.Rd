\name{gelman.plot}
\alias{gelman.plot}
\title{Gelman-Rubin-Brooks plot}

\usage{
gelman.plot(x, maxVars = 6, \dots)
} 

\arguments{
  \item{x}{A \code{posterior} object}
	\item{maxVars}{Maximum number of parameters to display.  The default
	  is small because this diagnostics can be time-consuming to calculate.}
  \item{...}{Additional arguments for \code{gelman.plot} in the \code{coda} package}
}

\description{
Creates a plot of Gelman and Rubin's shrink factor as a function of the number of iterations.
Computation is done by calling the \code{gelman.plot} function in the \code{coda} library.
See Brooks and Gelman (1998) or the help file for \code{gelman.plot} in the \code{coda} package
for details and more optional arguments.  The latter can be found in the Available Help section of the Help menu when 
	\code{coda} is loaded, which occurs when \code{gelman.plot} is called 
	on an object of class \code{posterior}.
}

\references{
Brooks, S. P. and A. Gelman (1997). General methods for monitoring
convergence of iterative simulations. \emph{Journal of Computational and
Graphical Statistics}, 7, 434-455.

Gelman, A. and D. B. Rubin (1992). Inference from iterative simulation
using multiple sequences. \emph{Statistical Science}, 7, 457-472.

}

\seealso{
   \code{\link{gelman.diag}},
   \code{\link{geweke.diag}},
   \code{\link{geweke.plot}}.
}
