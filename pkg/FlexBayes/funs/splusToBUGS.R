splusToBUGS <- function(data){

	if (is.list(data))
		dataBUGS <- data
	else
		dataBUGS <- list(data)
		
	len <- length(data)
	for (i in (1:len)){
		elem <- data[[i]]
		# If elem is a 3-D array, give a warning--not yet tested
		if (is.array(elem) && (length(dim(elem)) > 2))
			warning("3-D data arrays not yet tested in FlexBayes")
	  # If elem is a 2-D array, change from the data being stored
	  # by column to the data being stored by row.
		if (is.array(elem) && (length(dim(elem)) == 2)){
			elemVec <- as.vector(elem)
			elemBUGS <- matrix(data = elemVec, nrow = nrow(elem), 
				ncol = ncol(elem), byrow=T)
		} else {
			elemBUGS <- elem
		}
		dataBUGS[[i]] <- elemBUGS
	}
	return(dataBUGS)
}


