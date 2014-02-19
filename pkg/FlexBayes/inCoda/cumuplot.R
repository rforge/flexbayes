cumuplot.posterior  <- function(x, maxVars = 6, ...){
  if (nvar(x) > maxVars)
		x <- x[,1:maxVars]
	sims <- x$mcmc.list
	
	cumuplot(sims, ...)
}
