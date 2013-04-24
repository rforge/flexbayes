  ### FlexBayes prior distributions ###

bayes.normal <- function(mu = zero, sigma = identity, k0 = 1.0, ...)
{
  dots <- list(...)

  if(k0 <= 0)
    stop("parameter", sQuote("k0"), "is negative")

  if(!is.null(dots$mean.vector))
    mu <- dots$mean.vector

  if(!is.null(dots$covmat))
    sigma <- dots$covmat

  params <- list(mu = mu, sigma = sigma, k0 = k0)
  dstn <- list(name = "norm", label = "normal", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.t <- function(mu = zero, sigma = identity, df = 3, ...)
{
  dots <- list(...)

  if(!is.null(dots$mean.vector))
    mu <- dots$mean.vector

  if(!is.null(dots$covmat))
    sigma <- dots$covmat

	params <- list(mu = mu, sigma = sigma, df = df)
  dstn <- list(name = "t", label = "t", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.nonInformative <- function(mu = 0, sigma = 1, ...)
{
  dots <- list(...)

  if(!is.null(dots$mean.vector))
    mu <- dots$mean.vector

  if(!is.null(dots$covmat))
    sigma <- dots$covmat

	params <- list(mu = mu, sigma = sigma)
  dstn <- list(name = "nonInformative", label = "non-informative", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.beta <- function(shape1, shape2, ...)
{
	if(shape1 <= 0.0 || scale2 <= 0.0)
		stop("parameters must be positive")

	params <- list(shape1 = shape1, shape2 = shape2)
  dstn <- list(name = "beta", label = "beta", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.gamma <- function(shape, rate = 1.0, scale = 1.0/rate, ...)
{
	if(shape <= 0.0 || rate <= 0.0)
  	stop("parameters must be positive")

	params <- list(shape = shape, rate = rate)
  dstn <- list(name = "gamma", label = "gamma", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.normal.mixture <- function(mu = zero, sigma = identity, k = 3, props = c(0.5, 0.5), n.components = 2, k0 = 1.0, ...)
{
  dots <- list(...)

  if(!is.null(dots$mean.vector))
    mu <- dots$mean.vector

  if(!is.null(dots$covmat))
    sigma <- dots$covmat

  if(k0 <= 0)
    stop("parameter", sQuote("k0"), "is negative")

  if(k < 0)
    stop("parameter", sQuote("k"), "is negative")

  ## k = 0 indicates that covariances are different
  if(k == 0)
    n.components <- 2
       
	params <- list(mu = mu, sigma = sigma, k = k, props = props, components = n.components, k0 = k0)
  dstn <- list(name = "normmix", label = "normal mixture", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.t.mixture <- function(mu = zero, sigma = identity, k = 3, df = 3, props = c( 0.5, 0.5 ), n.components = 2, ...)
{
  dots <- list(...)

  if(!is.null(dots$mean.vector))
    mu <- dots$mean.vector

  if(!is.null(dots$covmat))
    sigma <- dots$covmat

  if(k < 0)
    stop("parameter", sQuote("k"), "is negative")

  ## k = 0 indicates that covariances are different
  if(k == 0)
    n.components <- 2
       
	params <- list(mu = mu, sigma = sigma, k = k, df = df, props = props, components = n.components )
  dstn <- list(name = "tmix", label = "t mixture", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.uniformShrinkage <- function(median = 1.0)
{
  if(median <= 0)
    stop(sQuote("median"), "is negative")

  params <- list(median = median)
  dstn <- list(name = "uniformShrinkage", label = "uniform shrinkage", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.duMouchel <- function(dispersion = 1.0)
{
  if(dispersion <= 0)
    stop(sQuote("dispersion"), "is negative")

  params <- list(dispersion = dispersion)
  dstn <- list(name = "duMouchel", label = "du Mouchel", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.nonInfoPower <- function(power = -1.0)
{
  params <- list(power = power)
  dstn <- list(name = "nonInfoPower", label = "non-informative power", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.massPoint <- function(value = 1.0)
{
  params <- list(value = value)
  dstn <- list(name = "massPoint", label = "mass point", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.wishart <- function(df = 1, scale = 1.0)
{
  if(df <= 0)
    stop("degrees of freedom is negative")

  params <- list(df = df, scale = scale)
  dstn <- list(name = "wishart", label = "Wishart", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.invWishart <- function(df = 1, scale = 1.0)
{
  if(df <= 0)
    stop("degrees of freedom is negative")

  params <- list(df = df, scale = scale)
  dstn <- list(name = "invWishart", label = "inverse Wishart", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


bayes.invChisq <- function( df = 3, sigma0.sq = 1)
{
	if(df < 0.0 )
  	stop("degrees of freedom is negative")

  if(sigma0.sq <= 0.0)
	  stop("scale is negative")
        
	params <- list(df = df, sigma0.sq = sigma0.sq)
  dstn <- list(name = "invChisq", label = "inverse chi-squared", parameters = params)
  oldClass(dstn) <- "fbdstn"
  dstn
}


print.fbdstn <- function(x, ...)
{
	cat(paste(x$label, "with:\n\n"))

	for(param in names(x$parameters)) {
		cat(paste(param, ":\n"))
		print(x$parameters[[param]])
		cat("\n")
	}

	invisible(x)
}



