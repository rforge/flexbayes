\name{raftery.diag}
\alias{raftery.diag}
\title{Raftery and Lewis's diagnostic}

\usage{raftery.diag(x, maxVars = 30, ...)}

\arguments{
  \item{x}{A \code{posterior} object}
	\item{maxVars}{Maximum number of parameters to display}
  \item{...}{Additional arguments for \code{raftery.diag} in the \code{coda} package}
}

\description{
  Raftery and Lewis's diagnostic.  Computation is done by calling the \code{raftery.diag} function in
  the \code{coda} library; see the \code{raftery.diag} section of the \code{coda} help
  for details and more optional arguments.  This can be found in the Available Help section of the Help menu when 
	\code{coda} is loaded, which occurs when \code{raftery.diag} is called 
	on an object of class \code{posterior}.
}
\seealso{
   \code{\link{geweke.diag}},
   \code{\link{gelman.diag}},
   \code{\link{heidel.diag}}.
}

