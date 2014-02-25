FlexBayesPriorDistributions <- c("normal", "t", "nonInformative", "beta",
                                 "gamma", "normal.mixture", "t.mixture",
                                 "uniformShrinkage", "duMouchel",
                                 "nonInfoPower", "massPoint", "wishart",
                                 "invWishart", "invChisq")

ParseDotsForParameters <- function(dots, params)
{
  dots.names <- names(dots)
  params.names <- names(params)

  if(is.null(dots.names)) {
    if(n.dots <- length(dots)) {
      n.dots <- min(c(length(params), n.dots))
      params[1:n.dots] <- dots
    }
  }

  else {
    named.dots <- dots.names[nchar(dots.names) > 0]
    if(!all(named.dots %in% params.names))
    stop("parameter not for dstn")

    common.names <- intersect(names(params), names(dots))
    params[common.names] <- dots[common.names]

    pidx <- !(params.names %in% common.names)
    didx <- !(dots.names %in% common.names)

    if((n.dots <- sum(didx)) > sum(pidx))
    stop("too many parameters")

    params[pidx][1:n.dots] <- dots[didx]
  }

  params
}


fbprior <- function(dstn, ...)
{
  dstn <- match.arg(dstn, choices = FlexBayesPriorDistributions)
  dots <- list(...)

  prior <- switch(dstn,

    "normal" = {
      default.parameters <- list(mean = 0, S = diag, k0 = 1)
      parameters <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = parameters)
    },

    "t" = {
      default.parameters <- list(mean = 0, S = diag, df = 1)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "nonInformative" = {
      default.parameters <- list(mean = 0, S = diag)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "beta" = {
      default.parameters <- list(shape1 = 0, shape2 = 0)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "gamma" = {
      default.parameters <- list(shape = 0, scale = 0)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "normal.mixture" = {
      default.parameters <- list(mean = 0, S = diag, k = 3, probs = rep(0.5, 2), k0 = 1)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "t.mixture" = {
      default.parameters <- list(mean = 0, S = diag, k = 3, df = 3, probs = rep(0.5, 2))
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "uniformShrinkage" = {
      default.parameters <- list(median = 1.0)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "duMouchel" = {
      default.parameters <- list(dispersion = 1.0)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "nonInfoPower" = {
      default.parameters <- list(power = -1.0)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "massPoint" = {
      default.parameters <- list(value = 1.0)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "wishart" = {
      default.parameters <- list(df = 1, Sigma = 1)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "invWishart" = {
      default.parameters <- list(df = 1, Sigma = 1)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "invChisq" = {
      default.parameters <- list(df = 3, sigma0.sq = 1)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    stop("Impossible error: switch in fbprior")
  )
  
  oldClass(prior) <- "fbprior"
  prior
}


print.fbprior <- function(x, ...)
{
	cat(paste(x$name, "with:\n\n"))

	for(name in names(x$parameters)) {
		cat(paste(name, ":\n"))
		print(x$parameters[[name]])
		cat("\n")
	}

	invisible(x)
}



