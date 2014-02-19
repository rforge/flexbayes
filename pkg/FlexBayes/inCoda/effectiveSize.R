# Calculate the effective sample sizes for a posterior object
effectiveSize.posterior  <- function(x, maxVars = 30){
  if (nvar(x) > maxVars)
		x <- x[,1:maxVars]
	sims <- x$mcmc.list

	effectiveSize(sims)
}


