#####################################################################################
##
##  MCMC CONTROL (internal)
##
#####################################################################################

sampler.control <- function(nBurnin, nSamples, nThin)
{
	bSize <- nBurnin
	simSize <- nSamples
	freqSize <- nThin
  
	if ( bSize <= 0 )
  {
	  warning( "The burn-in size is zero or negative, the default value of
    1000 was used instead." )
	  bSize <- 1000
	}
  
	if ( simSize <= 0 )
  {
	  warning( "The posterior sample size is zero or negative; the default
    value of 1000 was used instead" )
	  simSize <- 1000
	}
  
  if ( freqSize <= 0 )
  {
	  warning( "The posterior sample frequency is zero or negative; the default
    value of 1 was used instead" )
	  freqSize <- 1
	}
  
	list( bSize = bSize, simSize = simSize, freqSize = freqSize )
}


