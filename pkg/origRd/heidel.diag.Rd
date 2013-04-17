\name{heidel.diag}
\alias{heidel.diag}
\title{Heidelberger and Welch's diagnostic}

\usage{heidel.diag(x, maxVars = 30, ...)}

\arguments{
  \item{x}{A \code{posterior} object}
	\item{maxVars}{Maximum number of parameters to display}
  \item{...}{Additional arguments for \code{heidel.diag} in the \code{coda} package}
}

\description{
  Heidelberger and Welch's diagnostic.  Computation is done by calling the \code{heidel.diag}
  function of the \code{coda} library; see the \code{heidel.diag} section of the \code{coda} help
  for details and optional arguments.  This can be found in the Available Help section of the Help menu when 
	\code{coda} is loaded, which occurs when \code{heidel.diag} is called 
	on an object of class \code{posterior}.
}
\seealso{
   \code{\link{geweke.diag}},
   \code{\link{gelman.diag}},
   \code{\link{raftery.diag}}.
}
