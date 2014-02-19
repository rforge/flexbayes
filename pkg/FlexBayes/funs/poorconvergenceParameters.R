poorconvergenceParameters <- function (x, minESS = niter(x)*nchain(x) / 10){
  
	ess <- effectiveSize(x)
	poorConvergence1 <- (ess < minESS)
	
	if (nchain(x) > 1){
		# obtain the 97.5% quantile of the potential scale reduction factor
		gelmanRubin <- gelman.diag(x)$psrf[,2]
		poorConvergence2 <- (gelmanRubin > 1.2)
    
		# combine
		poorConvergenceParams <- varnames(x)[poorConvergence1 | poorConvergence2]
	} else {
	  poorConvergenceParams <- varnames(x)[poorConvergence1]
	}
	
	return (poorConvergenceParams)
}


