FlexBayesPriorDistributions <- c("normal", "t", "nonInformative", "beta",
                                 "gamma", "normalMixture", "tMixture",
                                 "uniformShrinkage", "duMouchel",
                                 "nonInfoPower", "massPoint", "wishart",
                                 "invWishart", "invChisq")

ParseDotsForParameters <- function(dots, params)
{
  dots.names <- names(dots)
  prp <- list()

  for(p in params) {
    if(is.null(prp[[p]] <- dots[[p]]))
      stop("parameter ", sQuote(p), " not found")
  }

  prp
}


fbprior <- function(dstn, ...)
{
  dstn <- match.arg(dstn, choices = FlexBayesPriorDistributions)
  dots <- list(...)

  prior <- switch(dstn,

    "normal" = {
      parameters <- ParseDotsForParameters(dots, c("mean", "S"))
      parameters$S <- as.matrix(parameters$S)
      if(!all(dim(parameters$S) == length(parameters$mean)))
        stop("mean vector ", sQuote("mean"), " and scale matrix ", sQuote("S"),
             " are not conformable")
      parameters$k0 <- 1.0
      list(name = dstn, parameters = parameters)
    },

    "t" = {
      parameters <- ParseDotsForParameters(dots, c("mean", "S", "df"))
      parameters$S <- as.matrix(parameters$S)
      if(!all(dim(parameters$S) == length(parameters$mean)))
        stop("mean vector ", sQuote("mean"), " and scale matrix ", sQuote("S"),
             " are not conformable")
      list(name = dstn, parameters = parameters)
    },

    "nonInformative" = {
      list(name = dstn, parameters = list())
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

    "normalMixture" = {
      parameters <- ParseDotsForParameters(dots, c("mean", "S", "w"))
      n.comps <- length(parameters$w)

      parameters$mean <- as.list(parameters$mean)
      p <- length(parameters$mean[[1]])
      if(any(sapply(parameters$mean, length) != p))
        stop("mean vectors are not all the same length")

      if(length(parameters$mean) != n.comps)
        stop("weights vector ", sQuote("w"), " and mean vectors ", sQuote("mean"),
             " are not conformable")

      parameters$S <- as.list(parameters$S)
      if(p > 1) {
        if(any(sapply(parameters$S, dim) != p))
          stop("mean vectors ", sQuote("mean"), " and scale matrix array ",
               sQuote("S"), " are not conformable")
      }

      else {
        if(any(sapply(parameters$S, length) != p))
           stop("mean vectors ", sQuote("mean"), " and scale matrix array ",
           sQuote("S"), " are not conformable")
      }

      if(length(parameters$S) != n.comps)
        stop("weights vector ", sQuote("w"), " and scale matrix array ",
             sQuote("S"), " are not conformable")

      parameters$k0 <- 1.0
      list(name = dstn, parameters = parameters)
    },

    "tMixture" = {
      parameters <- ParseDotsForParameters(dots, c("mean", "S", "w", "df"))
      n.comps <- length(parameters$w)

      parameters$mean <- as.list(parameters$mean)
      p <- length(parameters$mean[[1]])
      if(any(sapply(parameters$mean, length) != p))
        stop("mean vectors are not all the same length")

      if(length(parameters$mean) != n.comps)
        stop("weights vector ", sQuote("w"), " and mean vectors ", sQuote("mean"),
             " are not conformable")

      parameters$S <- as.list(parameters$S)
      if(p > 1) {
        if(any(sapply(parameters$S, dim) != p))
          stop("mean vectors ", sQuote("mean"), " and scale matrix array ",
               sQuote("S"), " are not conformable")
      }

      else {
        if(any(sapply(parameters$S, length) != p))
          stop("mean vectors ", sQuote("mean"), " and scale matrix array ",
               sQuote("S"), " are not conformable")
      }

      if(length(parameters$S) != n.comps)
        stop("weights vector ", sQuote("w"), " and scale matrix array ",
             sQuote("S"), " are not conformable")

      parameters$k <- 1.0
      list(name = dstn, parameters = parameters)
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
      default.parameters <- list(df = 1, scale = 1)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "invWishart" = {
      default.parameters <- list(df = 1, scale = 1)
      params <- ParseDotsForParameters(dots, default.parameters)
      list(name = dstn, parameters = params)
    },

    "invChisq" = {
      parameters <- ParseDotsForParameters(dots, c("df", "sigma0.sq"))
      list(name = dstn, parameters = parameters)
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
		cat(paste(name, ":\n", sep = ""))
		print(x$parameters[[name]])
		cat("\n")
	}

	invisible(x)
}



