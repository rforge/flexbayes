# Get the posterior samples for a particular chain
"[[.posterior" <- function (x, chainIndex){
	y <- posterior( sims = x$mcmc.list[[chainIndex]], DIC = x$DIC)
	y
}

# Subsetting function for posterior objects.  Parameters
# can be chosen by name using the post[,c("paramA", "paramB")] 
# syntax.
"[.posterior" <- function(x, i, j, drop=F){
  
	# If j is a list of parameter names,
	# check whether the names match existing parameters
	if (is.character(j)){
		paramNames <- dimnames( x$mcmc.list[[1]] )[[2]]
		if (!all(j %in% paramNames))
			stop("The parameters that are being selecting by name do not all exist in the posterior samples")
	}
	if (missing (i) && missing(j)){
		xSubset <- x
	} else if (missing(i)) {
		mcmcSubset <- list()
		nchains <- length( x$mcmc.list )
		for( k in (1:nchains) ){
		  thisMcmc <- x$mcmc.list[[k]]
		  attr( thisMcmc, "class" ) <- NULL
		  mcpar <- attr( thisMcmc, "mcpar" )
		  thisMcmc <- thisMcmc[,j,drop = F]
		  attr( thisMcmc, "class" ) <- "mcmc"
		  attr( thisMcmc, "mcpar" ) <- mcpar
		  mcmcSubset <- c( mcmcSubset, list( thisMcmc ) )
		}
		attr( mcmcSubset, "class" ) <- "mcmc.list"
		xSubset <- posterior (sims = mcmcSubset, DIC = x$DIC)
	} else if (missing(j)) {
		mcmcSubset <- list()
		nchains <- length( x$mcmc.list )

		for( k in (1:nchains) ){
		  thisMcmc <- x$mcmc.list[[k]]
		  attr( thisMcmc, "class" ) <- NULL
		  mcpar <- attr( thisMcmc, "mcpar" )
		  thisMcmc <- thisMcmc[i,,drop=F]
		  attr( thisMcmc, "class" ) <- "mcmc"
		  mcpar[1] <- 1  # set the starting iteration to 1
		  mcpar[2] <- length(i)  # set the ending iteration to length(i)
		  mcpar[3] <- 1  # set the thinning to 1
		  attr( thisMcmc, "mcpar" ) <- mcpar
		  mcmcSubset <- c( mcmcSubset, list( thisMcmc ) )
		}
		attr( mcmcSubset, "class" ) <- "mcmc.list"
		xSubset <- posterior (sims = mcmcSubset, DIC = x$DIC)
	} else {
		mcmcSubset <- list()
		nchains <- length( x$mcmc.list )

		for( k in (1:nchains) ){
		  thisMcmc <- x$mcmc.list[[k]]
		  attr( thisMcmc, "class" ) <- NULL
		  mcpar <- attr( thisMcmc, "mcpar" )
		  thisMcmc <- thisMcmc[i,j,drop=F]
		  attr( thisMcmc, "class" ) <- "mcmc"
		  mcpar[1] <- 1  # set the starting iteration to 1
		  mcpar[2] <- length(i)  # set the ending iteration to length(i)
		  mcpar[3] <- 1  # set the thinning to 1
		  attr( thisMcmc, "mcpar" ) <- mcpar
		  mcmcSubset <- c( mcmcSubset, list( thisMcmc ) )
		}
		attr( mcmcSubset, "class" ) <- "mcmc.list"
		xSubset <- posterior (sims = mcmcSubset, DIC = x$DIC)
  }
	return(xSubset)
}

# Window the samples in the posterior object--thin the samples,
# or choose a subset of the samples by specifying the indices of the new start and
# end iterations
"window.posterior" <- function(x, ...){
	loadCoda()
	windowed.mcmc.list <- window.mcmc.list(x$mcmc.list, ...)
	y <- posterior(sims = windowed.mcmc.list, DIC = x$DIC)
	y
}

