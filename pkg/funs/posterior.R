# Create a posterior object.  sims is the object of class mcmc.list,
# mcmc, matrix, data.frame that contains the posterior samples.  DIC 
# is the value of the DIC, estimated from the posterior samples.
posterior <- function(sims, DIC = NULL, call = NULL){

	if (!(class(sims) %in% c("mcmc.list", "mcmc", "matrix", "data.frame")))
		stop ("Creating an object of class posterior requires an object of class mcmc, mcmc.list, matrix, or data.frame containing the posterior samples")
	if (class(sims) != "mcmc.list")
		simsList <- mcmc.list(mcmc(sims))
	else {
	  # If there is only one parameter and the mcmc objects are vectors, change them to matrices
	  mcpar <- attr(sims[[1]], "mcpar")
	  for (i in (1:length(sims))){
		  if( !is.matrix(sims[[i]]) ){
		    sims[[i]] <- matrix(as.vector(sims[[i]]), ncol=1)
		    dimnames(sims[[i]]) <- list(NULL, c(""))
		  }
		  attr(sims[[i]], "mcpar") <- mcpar
		  attr(sims[[i]], "class") <- "mcmc"
	  }

		simsList <- sims
  }
	post <- list (mcmc.list = simsList, DIC = DIC, call = call) 
	attr(post, "class") <- "posterior"	
	post
}
