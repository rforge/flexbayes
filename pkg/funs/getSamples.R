getSamples <- function(x, param = NULL, chainIndex = NULL){
  if( is.character( param ) && 
    !all( is.element( param, varnames( x ) ) ) )
    stop("Some of the parameters selected by name do not exist in the posterior object")
  if (is.null(chainIndex)){
  	firstChain <- x$mcmc.list[[1]]
  	## change the class of firstChain,
  	## which is currently an mcmc object.
    attr( firstChain, "class" ) <- NULL
  	if (is.null(param))
      samples <- firstChain
    else {
  	  samples <- firstChain[,param]
  	}
  	nchains <- nchain(x)
  	if (nchains > 1)
  		for (i in 2:nchains){
  			thisChain <- x$mcmc.list[[i]]
  			if (is.null(param)){
  			  # concatenate the columns of thisChain 
  			  # onto the columns of samples
  			  if( nvar( x ) == 1 )
  			    samples <- c(samples, thisChain)
  			  else
    			  samples <- rbind(samples, thisChain)
  			} else {
  			  # change the class of thisChain so that
  			  # the following subset operator will not
  			  # call the mcmc subset function if coda
  			  # is loaded
  			  attr( thisChain, "class" ) <- NULL
  			  if( length( param ) == 1 )
  			    samples <- c(samples, thisChain[,param])
  			  else
    			  samples <- rbind(samples, thisChain[,param])
  			}
	  	}
  } else {
    if (is.null(param)){
      samples <- as.matrix( x$mcmc.list[[chainIndex]] )
      if( nvar( x ) == 1 )
        samples <- as.vector( samples )  
    } else {
      samples <- as.matrix( x$mcmc.list[[chainIndex]][,param] )
      if( length(param) == 1 )
        samples <- as.vector( samples )  
    }
  }
  return(samples)
}
