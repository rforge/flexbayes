\name{crosscorr.posterior}
\alias{crosscorr}
\alias{crosscorr.posterior}
\title{Calculates the correlation matrix for a posterior object}

\usage{crosscorr(x, maxVars = 30)}

\arguments{
\item{x}{A \code{posterior} object}
\item{maxVars}{Maximum number of parameters to display}
}

\description{
	Calculates the pairwise correlation matrix for the parameters in a \code{posterior} object.  
	If there are multiple chains then the samples from all of the chains are pooled for the calculation.
	Computation is done by calling the \code{crosscorr} function in the \code{coda} library.  
}
\value{
  A matrix containing the pairwise correlations of the parameters.
}
\seealso{
\code{\link{autocorr}},
\code{\link{lagged.crosscorr}}.
}
