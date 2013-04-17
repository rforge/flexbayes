# Geweke diagnostic
geweke.diag.posterior  <- function(x, maxVars = 30, ...){
  if (nvar(x) > maxVars)
		x <- x[,1:maxVars]
	sims <- x$mcmc.list

	geweke.diag(sims, ...)
}

geweke.plot.posterior  <- function(x, maxVars = 6, ...){
  if (nvar(x) > maxVars)
		x <- x[,1:maxVars]
	sims <- x$mcmc.list

	geweke.plot(sims, ...)
}
