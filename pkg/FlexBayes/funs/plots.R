# This function is intended to create plots summarizing the posterior distributions
# of the parameters.  However, it is not yet working.

#plotPosterior <- function(x, ...){
#	sims <- x$mcmc.list
#	nChains <- length(sims)
#	nSimsPerChain <- dim(sims[[1]])[1]
#  nParams <- dim(sims[[1]])[2]
#
#	paramNames <- dimnames(sims[[1]])[[2]]
#	simsArray <- array(data=0, dim = c(nSimsPerChain, nChains, nParams),
#		dimnames = list(c(),c(),rev(paramNames)))
#	for (i in (1:nChains))
#		for (j in (1:nParams))
#			simsArray[,i,nParams-j+1] <- sims[[i]][,j]
#
#	bugsArray <- as.bugs.array (simsArray, DIC = x$is.DIC)
#
#	plot(bugsArray)
#}

# Create traceplots for a posterior object
traceplot.posterior <- function(x, maxVars = 30, ...){
	
	if (nvar(x) > maxVars)
		x <- x[,1:maxVars]
	sims <- x$mcmc.list
	
	traceplot(sims, maxVars, ...)
}

# Create autocorrelation plots for a posterior object
autocorr.plot.posterior <- function(x, maxVars = 30, ask = F, ...){
		
	if (nvar(x) > maxVars)
		x <- x[,1:maxVars]
	sims <- x$mcmc.list

	invisible(autocorr.plot(x=sims, maxVars=maxVars, ask=ask, ...))

}

# Creates a cross-correlation plot for a posterior object
crosscorr.plot.posterior <- function(x, maxVars = 30, ...){
		
	if (nvar(x) > maxVars)
		x <- x[,1:maxVars]
	sims <- x$mcmc.list
	
	crosscorr.plot(x=sims, maxVars=maxVars, ...)

}

