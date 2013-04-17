# Calculate the lagged autocorrelations by parameter for a posterior object
autocorr.posterior <- function (x, maxVars = 30, ...) 
{
  if (nvar(x) > maxVars)
		x <- x[,1:maxVars]
	sims <- x$mcmc.list
	
	loadCoda()
	autocorr.diag(mcmc.obj = sims, ...)
}

# Calculate the lagged crosscorrelations for a posterior object
lagged.crosscorr.posterior <- function (x, maxVars = 30, ...) 
{
	if (nvar(x) > maxVars)
		x <- x[,1:maxVars]
	sims <- x$mcmc.list
	
	lagged.crosscorr(x=sims, ...)
}

# Calculate the parameter cross-correlations for a posterior 
# object
crosscorr.posterior <- function(x, maxVars = 30){
	if (nvar(x) > maxVars)
		x <- x[,1:maxVars]
	sims <- x$mcmc.list
	
	crosscorr(sims)
}

