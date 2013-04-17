\name{summary.posterior}
\alias{summary}
\alias{summary.posterior}
\title{Summarize posterior inferences}

\usage{summary(x, maxVars=30, digits=4)}

\arguments{
\item{x}{A \code{posterior} object}
\item{maxVars}{Maximum number of parameters to display}
\item{digits}{The number of significant digits to display}
}

\description{
	Obtains posterior means, standard deviations, relevant quantiles,
	and DIC if available.
}
\seealso{
   \code{\link{nvar.posterior}},
   \code{\link{start.posterior}},
   \code{\link{varnames.posterior}},
   \code{\link{subset.posterior}}.
}

