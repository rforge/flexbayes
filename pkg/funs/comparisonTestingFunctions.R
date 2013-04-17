compareSampleDistributions <- function (samples1, samples2, 
	print.ks=T, pvalue.cutoff=0.01){
#######################################################################
##
## compareSampleDistributions
##
## compare the distributions of two vectors of samples using a 
## Kolmogorov-Smirinov test.  If print.ks = T, print the output of the K-S test.
## return T if the pvalue is greater than pvalue.cutoff (sample 
## distributions are indistinguishable).  Return an indicator of whether
## the pvalue for the test was greater than pvalue.cutoff (whether the 
## test passed).
## 

	nvars <- ncol(samples1)
	if (nvars != ncol(samples2))
		stop("the number of variables should be equal in the two data sets")
	ks.test <- list(-1, nvars)
	p.value <- rep(-1, nvars)
	
	for (i in (1:nvars)){
		samples1.thisVar <- samples1[,i]
		samples2.thisVar <- samples2[,i]
		ks.test[[i]] <- ks.gof(samples1.thisVar, samples2.thisVar)
		p.value[i] <- ks.test[[i]]$p.value
	}
	
	if (print.ks)
		print( ks.test )
		
	return ( all( p.value >= ( pvalue.cutoff / nvars ) ) )	
}
