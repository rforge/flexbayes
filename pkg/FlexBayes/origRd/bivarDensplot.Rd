\name{bivarDensplot}
\alias{bivarDensplot}
\title{Plot bivariate posterior densities}

\usage{
bivarDensplot( post, maxVars = 5, cuts = 8,
  xylimProb = 0.98, lwd = 1.5, col = 1, ...)
}
\arguments{
	\item{post}{An object of class \code{posterior}}
	\item{maxVars}{Maximum number of parameters to display}
	\item{cuts}{The number of contour lines on the plot}
	\item{xylimProb}{Bounds on the display region, in terms of the posterior probability.  
	  For each parameter, the axis limits are equal to the bounds of a posterior interval 
	  with probability xylimProb.  Decreasing this value zooms in on the highest-probability 
	  region.}
	\item{lwd}{The line width on the displayed plots.} 
	\item{col}{The color of lines on the displayed plots.}
	\item{...}{Arguments for the \code{strip.default} function, which is called when drawing 
	  the individual plot titles.}
}
\description{
	Creates contour plots of the estimated bivariate posterior densities of the parameters.  The samples
	in all of the chains are pooled when obtaining these estimates.  
}
\details{
  Due to the smoothing inherent in kernel density estimates, the density can be estimated to be positive in regions that are not in the 
  support of the posterior density.  For instance, the kernel density estimate of a variance parameter can be positive just below zero.  
  Also, in the \code{graphsheet} device the individual plot titles can sometimes appear too low, so that they run into the line below.
  To avoid this issue the user can open the \code{java.graph} device instead, by calling \code{java.graph()} before calling 
  the \code{bivarDensplot} function.
}

\seealso{
   \code{\link{summary.posterior}}.
}

\examples{
##########################################
# Fit an overdispersed Poisson model to the
# pumps data, then plot the resulting 
# bivariate density estimates

xi.prior = bayes.uniformShrinkage (0.5)
random.var.prior = bayes.invWishart( df = 3, scale = diag( c(1,1) ) )

pump.prior = bhpm.prior ( xi = xi.prior,
  random.var =  random.var.prior, 
  common.glm = 2 )

# Specify the control parameters (burn-in length, 
# etc.) and the initial values for the parameters.  
pump.sampler <- bhpm.sampler( nSamples=2000, nThin= 100, 
  nChains = 1, nBurnin = 1000 ) 

pump.exposure <- ~ e
pump.random <- z ~ 1 + x

## call bhpm to fit a gamma-conjugate model
pump.bhpm <- bhpm( random.formula = pump.random, 
  exposure.formula = pump.exposure, data = pumps, 
  prior = pump.prior, sampler = pump.sampler,
  overdispersion = "gamma-conj" )

# check convergence
autocorr.plot(pump.bhpm)

# plot the bivariate density estimates
bivarDensplot(pump.bhpm[,6:10])

}

