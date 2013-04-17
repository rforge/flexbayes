# Calculate the highest posterior density intervals for a posterior object
HPDinterval.posterior  <- function(x, maxVars = 30,...){
  if (nvar(x) > maxVars)
		x <- x[,1:maxVars]
	sims <- x$mcmc.list

	HPDinterval(sims, ...)
}
