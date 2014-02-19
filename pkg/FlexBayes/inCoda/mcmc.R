# functions to create mcmc and mcmc.list objects.

mcmc <- function(data, start=1, end=-1, thin=1, burnin=0){
	if (is.data.frame(data))
		data <- as.matrix(data)
	if (is.vector(data))
		data <- matrix(data, ncol = 1)
	if ( !(is.matrix(data)) )
		stop("data must be a matrix, data frame, or vector.")
	if (!all(is.numeric(as.vector(data))))
		stop("data must be an all-numeric matrix")
	
	if (start == 1)
    start <- start + burnin
	if (end == -1){
		end <- ( dim(data)[1] - 1 ) * thin
		end <- end + start
  }
	
	attr(data, "class") <- "mcmc"
	attr(data, "mcpar") <- c(start, end, thin)

	data
}

is.mcmc.list <- function(x){
	return (!is.null(attr(x, "class")) && (attr(x,"class") == "mcmc.list"))
}

is.mcmc <- function(x){
	return (!is.null(attr(x, "class")) && (attr(x,"class") == "mcmc"))
}

mcmc.list <- function(...){
	mcmcList <- list(...)
	if( is.list( mcmcList[[1]] ) )
	  mcmcList <- mcmcList[[1]]

  for( i in (1:length(mcmcList)) ){
    if( !is.mcmc( mcmcList[[i]] ) )
      stop("Cannot create an mcmc.list object from non-mcmc components");
  }
	attr(mcmcList, "class") <- "mcmc.list"
	return(mcmcList)
}
