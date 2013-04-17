# Heidelberger and Welch's diagnostic
heidel.diag.posterior  <- function(x, maxVars = 30, ...){
  if (nvar(x) > maxVars)
		x <- x[,1:maxVars]
	sims <- x$mcmc.list

	heidel.diag(sims, ...)
}
