# change the functions from CODA so that calls to the function
# will route to the ".posterior" version if the object is of class posterior,
# and to the original version of the object is of class mcmc or mcmc.list

"crosscorr" <- function(x, maxVars = 30){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		crosscorr.mcmc <- get(name = "crosscorr", where = "coda")
		crosscorr.mcmc(x)
	} else if (class(x) == "posterior")
		crosscorr.posterior(x, maxVars)
}

"effectiveSize" <- function(x, maxVars = 30){
	
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		
		loadCoda()	 
		effectiveSize.mcmc <- get(name = "effectiveSize", where = "coda")
		
		# the call to effectiveSize.mcmc creates unnecessary warnings,
		# due to a problem with the ar function.
		# turn off these warnings.
		oldWarn <- options("warn")$warn
	  options(warn=-1)
	  on.exit(options(warn=oldWarn))
	  
	  round(effectiveSize.mcmc(x))
	  
	} else if (class(x) == "posterior")
		round(effectiveSize.posterior(x, maxVars))
}

"autocorr.plot" <- function(x, maxVars = 30, ...){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		autocorr.plot.mcmc <- get(name = "autocorr.plot", where = "coda")
		invisible(autocorr.plot.mcmc(x, ...))
	} else if (class(x) == "posterior"){
		invisible(autocorr.plot.posterior(x, maxVars, ...))
	}
}

"crosscorr.plot" <- function(x, maxVars = 30, ...){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		crosscorr.plot.mcmc <- get(name = "crosscorr.plot", where = "coda")
		invisible(crosscorr.plot.mcmc(x, ...))
	} else if (class(x) == "posterior"){
		invisible(crosscorr.plot.posterior(x, maxVars, ...))
	}
}

"traceplot" <- function(x, maxVars = 30, ...){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		traceplot.mcmc <- get(name = "traceplot", where = "coda")
		traceplot.mcmc(x, ...)
	} else if (class(x) == "posterior"){
		traceplot.posterior(x, maxVars, ...)
	}
}

"lagged.crosscorr" <- function(x, maxVars = 30, ...){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		lagged.crosscorr.mcmc <- get(name = "autocorr", where = "coda")
		lagged.crosscorr.mcmc(x, ...)
	} else if (class(x) == "posterior"){
		lagged.crosscorr.posterior(x, maxVars, ...)
	}
}

"autocorr" <- function(x, maxVars = 30, ...){
	if ((class(x) == "mcmc.list")||(class(x) == "mcmc")){
		loadCoda()	 
		autocorr.mcmc <- get(name = "autocorr", where = "coda")
		autocorr.mcmc(x, ...)
	} else if (class(x) == "posterior"){
		autocorr.posterior(x, maxVars, ...)
	}
}

"densplot" <- function(x, maxVars = 30, ...){
  if ((class(x) == "mcmc.list")||(class(x) == "mcmc")){
		loadCoda()	 
		densplot.mcmc <- get(name = "densplot", where = "coda")
		densplot.mcmc(x, ...)
	} else if (class(x) == "posterior"){
		invisible(densplot.posterior(x, maxVars = maxVars, ...))
	}
}
	
"nvar" <- function(x){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()
		nvar.mcmc <- get(name = "nvar", where = "coda")
		nvar.mcmc(x)
	} else if (class(x) == "posterior")
		nvar.posterior(x)
		
}

"thin" <- function(x){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()
		thin.mcmc <- get(name = "thin", where = "coda")
		thin.mcmc(x)
	} else if (class(x) == "posterior")
		thin.posterior(x)
}

"varnames" <- function(x, allow.null=T){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		varnames.mcmc <- get(name = "varnames", where = "coda")
		varnames.mcmc(x, allow.null)
	} else if (class(x) == "posterior")
		varnames.posterior(x)
}
	
"varnames<-" <- function(x, value){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		varnames.mcmc <- get(name = "varnames<-", where = "coda")
		varnames.mcmc(x, value)
	} else if (class(x) == "posterior"){
		varnames.posterior(x) <- value
		x
  }
}

"nchain" <- function(x){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		nchain.mcmc <- get(name = "nchain", where = "coda")
		nchain.mcmc(x)
	} else if (class(x) == "posterior")
		nchain.posterior(x)
}

"niter" <- function(x){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		niter.mcmc <- get(name = "niter", where = "coda")
		niter.mcmc(x)
	} else if (class(x) == "posterior")
		niter.posterior(x)
}

"as.mcmc.list" <- function(x){
	if (class(x) == "posterior")
		as.mcmc.list.posterior(x)
  else {
		loadCoda()	 
		as.mcmc.list.mcmc <- get(name = "as.mcmc.list", where = "coda")
		as.mcmc.list.mcmc(x)
	}
}

"gelman.diag" <- function(x, maxVars = 30, ...){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		gelman.diag.mcmc <- get(name = "gelman.diag", where = "coda")
		gelman.diag.mcmc(x, ...)
	} else if (class(x) == "posterior")
		gelman.diag.posterior(x, maxVars, ...)
}

"gelman.plot" <- function(x, maxVars = 6, ...){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		gelman.plot.mcmc <- get(name = "gelman.plot", where = "coda")
		invisible(gelman.plot.mcmc(x, ...))
	} else if (class(x) == "posterior")
		invisible(gelman.plot.posterior(x, maxVars, ...))
}

"geweke.diag" <- function(x, maxVars = 30, ...){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		geweke.diag.mcmc <- get(name = "geweke.diag", where = "coda")
		geweke.diag.mcmc(x, ...)
	} else if (class(x) == "posterior")
		geweke.diag.posterior(x, maxVars, ...)
}

"geweke.plot" <- function(x, maxVars = 6, ...){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		geweke.plot.mcmc <- get(name = "geweke.plot", where = "coda")
		invisible(geweke.plot.mcmc(x, ...))
	} else if (class(x) == "posterior")
		invisible(geweke.plot.posterior(x, maxVars, ...))
}

"raftery.diag" <- function(x, maxVars = 30, ...){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		raftery.diag.mcmc <- get(name = "raftery.diag", where = "coda")
		raftery.diag.mcmc(x, ...)
	} else if (class(x) == "posterior")
		raftery.diag.posterior(x, maxVars, ...)
}

"heidel.diag" <- function(x, maxVars = 30, ...){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		heidel.diag.mcmc <- get(name = "heidel.diag", where = "coda")
		heidel.diag.mcmc(x, ...)
	} else if (class(x) == "posterior")
		heidel.diag.posterior(x, maxVars, ...)
}

"cumuplot" <- function(x, maxVars = 6, ...){
	if ((class(x) == "mcmc.list") || (class(x) == "mcmc")){
		loadCoda()	 
		cumuplot.mcmc <- get(name = "cumuplot", where = "coda")
		cumuplot.mcmc(x, ...)
	} else if (class(x) == "posterior")
		cumuplot.posterior(x, maxVars, ...)
}
