##
#
#	Distributions for bayes priors
#
##

bayes.normal <- function( mean.vector = zero, covmat = identity, k0 = 1 )
{

        if ( k0 <= 0 )
        {
          stop( "bayes.normal: Parameter k0 must be positive." )
        }

	name <- "normal"
	params <- list( "mean vector" = mean.vector, "covariance matrix" = covmat, "k0" = k0 )

	new("bayes.distribution", name = name, parameters = params)
}


bayes.t <- function(mean.vector = zero, covmat = identity, df = 3 )
{
	name <- "t"
	params <- list("mean vector" = mean.vector, "covariance matrix" = covmat,
		"t degrees of freedom" = df)

	new("bayes.distribution", name = name, parameters = params)
}


##added for compatibility
bayes.nonInformative <- function( mean.vector = 0, covmat = 1 )
{
	name <- "non-informative"
	params <- list( "mean vector" = mean.vector, "covariance matrix" = covmat )

	new( "bayes.distribution", name = name, parameters = params )
}


bayes.beta <- function(shape = 0, scale = 0 )
{
	if(shape <= 0.0 || scale <= 0.0)
		stop("the shape and scale must both be positive.")

	name <- "beta"
	params <- list(shape = shape, scale = scale)

	new("bayes.distribution", name = name, parameters = params)
}


bayes.gamma <- function(shape = 0, scale = 0)
{
	if(shape <= 0.0 || scale <= 0.0)
		stop("the shape and scale must both be positive.")

	name <- "gamma"
	params <- list(shape = shape, scale = scale)

	new("bayes.distribution", name = name, parameters = params)
}


bayes.normal.mixture <- function( mean.vector = zero, covmat = identity, k = 3, props = c( 0.5, 0.5 ), n.components = 2, k0 = 1  )
{
        ## k = 0 indicates that covariances are different

        if ( k < 0 )
        {
          stop( "bayes.normal.mixture: Parameter k must be positive." )
        }
        else if ( k == 0 )
        {
          n.components <- 2
        }
       
	name <- "normal mixture"
	params <- list( "mean vector" = mean.vector, "covariance matrix" = covmat, "k" = k, "props" = props, "components" = n.components, "k0" = k0 )

	new("bayes.distribution", name = name, parameters = params)
}


bayes.t.mixture <- function( mean.vector = zero, covmat = identity, k = 3, df = 3, props = c( 0.5, 0.5 ), n.components = 2 )
{
        ## k = 0 indicates that covariances are different

        if ( k < 0 )
        {
          stop( "bayes.t.mixture: Parameter k must be positive." )
        }
        else if ( k == 0 )
        {
          n.components <- 2
        }
       
	name <- "t mixture"
	params <- list( "mean vector" = mean.vector, "covariance matrix" = covmat, "k" = k, "t degrees of freedom" = df, "props" = props, "components" = n.components )

	new("bayes.distribution", name = name, parameters = params)
}



bayes.uniformShrinkage <- function( median.value = 1.0 )
{
  if ( median.value <= 0 )
  {
    stop( "bayes.uniformShrinkage: median.value must be a positive number." )
  }

  name <- "uniform shrinkage"
  params <- list( "median" = median.value  )

  new( "bayes.distribution", name = name, parameters = params )

}#end



bayes.duMouchel <- function( dispersion = 1.0 )
{
  if ( dispersion <= 0 )
  {
    stop( "bayes.uniformShrinkage: dispersion must be a positive number." )
  }

  name <- "du Mouchel"
  params <- list( "dispersion" = dispersion )

  new( "bayes.distribution", name = name, parameters = params )

}#end



bayes.nonInfoPower <- function( power = -1.0 )
{
  name <- "non-informative power"
  params <- list( "power" = power )

  new( "bayes.distribution", name = name, parameters = params )

}#end



bayes.massPoint <- function( value = 1 )
{
  name <- "mass point"
  params <- list( "value" = value )

  new( "bayes.distribution", name = name, parameters = params )

}#end



bayes.wishart <- function( df = 1, scale = 1 )
{

  if ( df <= 0 )
  {
    stop( "bayes.wishart: degrees of freedom must be positive." )
  }

  name <- "wishart"
  params <- list( "degrees of freedom" = df, "scale" = scale )

  new( "bayes.distribution", name = name, parameters = params )

}#end
  

bayes.invWishart <- function( df = 1, scale = 1 )
{

  if ( df <= 0 )
  {
    stop( "bayes.invWishart: degrees of freedom must be positive." )
  }

  name <- "invWishart"
  params <- list( "degrees of freedom" = df, "scale" = scale )

  new( "bayes.distribution", name = name, parameters = params )

}#end
  
bayes.invChisq <- function( df = 3, sigma0.sq = 1)
{
	if ( any( df < 0.0 ) || any( sigma0.sq <= 0.0 ) ){
	  stop( "the df and scale vectors must both be positive." )
  }
        
	name <- "invChisq"
	params <- list( df = df, sigma0.sq = sigma0.sq )

	new( "bayes.distribution", name = name, parameters = params )
}

print.bayes.distribution <- function(object, digits)
{
	name <- object@name
	params <- object@parameters

	cat(paste(name, "with:\n\n"))

	for(name in names(params)) {
		cat(paste(name, ":\n"))
		print(params[[name]])
		cat("\n")
	}
	invisible(object)
}



