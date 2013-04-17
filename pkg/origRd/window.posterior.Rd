\name{window.posterior}
\alias{window.posterior}
\alias{window}
\title{Windowing for a posterior object}

\usage{window(x, start = start(x), end = end(x), thin = 1)}

\arguments{
\item{x}{A \code{posterior} object}

\item{start}{The first desired iteration.  Should be no smaller than \code{start(x)}.
	If it is smaller, then \code{start} will be set to \code{start(x)}.}	

\item{end}{The last desired iteration.  Should be no larger than \code{end(x)}.
	If it is larger, then \code{end} will be set to \code{end(x)}.}	

\item{thin}{The desired thinning (lag between saved parameter values)}
}

\description{
	Windowing function for \code{posterior} objects.  The samples can be thinned by the 
	desired amount, or a subset of the iterations of the chains can be chosen using \code{start}
	and \code{end}.
	}

