
# Get the names of the parameters in a posterior object
varnames.posterior <- function(x){
  if( is.null( dimnames(x$mcmc.list[[1]]) ) )
    varnameVec <- NULL
  else
    varnameVec <- dimnames(x$mcmc.list[[1]])[[2]]
	return(varnameVec)
}

# Set the names of the parameters in a posterior object
"varnames.posterior<-" <-
function (x, value) {
  temp <- dimnames(x$mcmc.list[[1]])
  temp[[2]] <- value
  dimnames(x$mcmc.list[[1]]) <- temp
  x
}

# Get the number of parameters in a posterior object
nvar.posterior <- function(x){
	return( dim( x$mcmc.list[[1]] )[2] )
}

# Get the number of chains in a posterior object
nchain.posterior <- function(x){
	return(length(x$mcmc.list))
}

# Get the number of iterations in a posterior object
niter.posterior <- function(x){
	return(dim(x$mcmc.list[[1]])[1])
}

# Get the index of the first iteration that was saved in
# the MCMC run when creating x (a posterior object)
"start.posterior" <- function(x){
	firstMcmc <- x$mcmc.list[[1]]
	return(attr(firstMcmc, "mcpar")[1])
}

# Get the index of the last iteration that was saved in
# the MCMC run when creating x (a posterior object)
"end.posterior" <- function(x){
	firstMcmc <- x$mcmc.list[[1]]
	return(attr(firstMcmc, "mcpar")[2])
}

# Get the amount of thinning in
# the MCMC run when creating x (a posterior object)
"thin.posterior" <- function(x){
	firstMcmc <- x$mcmc.list[[1]]
	return(attr(firstMcmc, "mcpar")[3])
}

# Get the parameter samples for chain chainIndex as an mcmc.list object
as.mcmc.list.posterior <- function(x){
	x$mcmc.list
}
