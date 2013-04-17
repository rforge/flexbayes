\name{densplot.posterior}
\alias{densplot}
\alias{densplot.posterior}
\title{Plot parameter posterior densities}

\usage{
densplot(x, maxVars = 30, maxPerPage = 6, 
  region = "credible", kernel.width = "nrd", 
  level = 0.95, lwd = 2, col = 1, 
  cex.axis = 0.5, cex = 0.7, ...)
}
\arguments{
	\item{x}{An object of class \code{posterior}}
	\item{maxVars}{Maximum number of parameters to display}
	\item{maxPerPage}{Maximum number of plots per page.  Can be 1 or any
	  even number (the number of columns is fixed at 2 for maxPerPage > 1).}
	\item{region}{The type of posterior interval to shade on the plots.  One
	  of \code{"hpd"}, \code{"credible"}, or \code{"none"}, yielding a highest posterior density region,
	  a credible interval such that the probability mass above the interval is the same
	  as below, or no shaded region, respectively.}
	\item{kernel.width}{Can be either a numeric value, specifying the width of the kernel in the 
	  density estimation, or a string specifying the method for computing this width.  See
	  the documentation of the \code{width} argument to the \code{density} function for more 
	  information.}
	\item{level}{The probability coverage of the credible interval that will be displayed.}
	  E.g., \code{level = 0.95} specifies that a 95% credible interval should be displayed.}
	\item{lwd}{The line width on the displayed plots.} 
	\item{col}{The color of lines on the displayed plots.}
	\item{cex.axis}{Character expansion of tick mark labels relative to device's standard size.}
	\item{cex}{Character expansion of plot titles relative to device's standard size.}
	\item{...}{Arguments for the \code{strip.default} function, which is called when drawing 
	  the individual plot titles.}
}
\description{
	Plots the estimated posterior densities of the parameters.  For each parameter the samples
	in all of the chains are pooled in order to obtain a single estimated density function.  
	The plots can have shaded posterior intervals.
}
\details{
  Due to the smoothing inherent in kernel density estimates, the density can be estimated to be positive in regions that are not in the 
  support of the posterior density.  For instance, the kernel density estimate of a variance parameter can be positive just below zero.
}

\seealso{
   \code{\link{summary.posterior}}.
}

\examples{
##########################################
# Fit an overdispersed Poisson model to the
# pumps data, then plot the resulting 
# posterior density estimates

xi.prior = bayes.uniformShrinkage (0.5)
random.var.prior = bayes.invWishart( df = 3, scale = diag( c(1,1) ) )

pump.prior = bhpm.prior ( xi = xi.prior,
  random.var =  random.var.prior, 
  common.glm = 2 )

# Specify the control parameters (burn-in length, 
# etc.) and the initial values for the parameters.  
pump.sampler <- bhpm.sampler( nSamples=1000, nThin= 50, 
  nChains = 1, nBurnin = 1000,
  init.point = "prior", update.cov = 0 ) 

pump.exposure <- ~ e
pump.random <- z ~ 1 + x

## call bhpm to fit a gamma-conjugate model
pump.bhpm <- bhpm( random.formula = pump.random, 
  exposure.formula = pump.exposure, data = pumps, 
  prior = pump.prior, sampler = pump.sampler,
  overdispersion = "gamma-conj" )

# check convergence
autocorr.plot(pump.bhpm)

# plot the posterior density estimates
densplot(pump.bhpm)

}

