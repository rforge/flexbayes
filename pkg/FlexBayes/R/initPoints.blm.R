# return: draws from beta and sigma2 based on mixture of prior and MLE's
#         (actual draws are returned by generateInitPoints.blm() below).
#         (useful for creating different chains starting at the draws for the
#         Bayes Linear Model)

init.points.blm <- function(formula, data, prior = blm.prior(),
                            number.draws = 1, na.action = na.fail,
                            contrasts = NULL, ...)
{
	model.blm <- call("model.frame", formula = formula, na.action = na.action)

	if(!missing( data))
		model.blm$data <- data

  model.blm <- eval(model.blm, parent.frame())
	Terms <- attr(model.blm, "terms")
	Y <- model.extract(model.blm, "response")
	X <- model.matrix(Terms, model.blm, contrasts)

	dim.Cov <- dim(X)[2]
	number.data <- length(Y)

	## Validate/interpret the arguments specifying the prior ##
	if(class(prior$priorBeta) == "fbdstn" )
	  prior.beta <- prior$priorBeta$name

  else
	  prior.beta <- "non-informative"

	if(class(prior$priorSigma) == "fbdstn")
	  prior.sigma <- prior$priorSigma$name

	else
	  prior.sigma <- "non-informative"

  beta.components <- 1
  beta.props <- 1

	iType <- paste(prior.beta, prior.sigma, sep = ":")

	iType <- switch(iType,

		"non-informative:non-informative" = {
			betaMean <- rep(0.0, dim.Cov)
			betaCov <- diag( dim.Cov )
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			0
		},

		"non-informative:invChisq" = {
			betaMean <- rep(0.0, dim.Cov)
			betaCov <- diag(dim.Cov)
			betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			0
		},

		"normal:non-informative" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			1
		},

		"normal:invChisq" = {
			betaMean <- prior$priorBeta@parameters[["mu"]]
			betaCov <- prior$priorBeta@parameters[["sigma"]]
			betaDf <- 3.0 
			sigmaDf <- prior$priorSigma@parameters[["df"]]
			sigmaScale <- prior$priorSigma@parameters[["sigma0.sq"]]
			ifelse(prior$conjugate, 6, 1)
		},


		"normmix:non-informative" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
      beta.k <- prior$priorBeta@parameters[["k"]]
      beta.props <- prior$priorBeta@parameters[["props"]]
      beta.components <- prior$priorBeta@parameters[["components"]]

			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			7
		},

    "normmix:invChisq" = {
      betaMean <- prior$priorBeta$parameters[["mu"]]
      betaCov <- prior$priorBeta$parameters[["sigma"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["props"]]
      beta.components <- prior$priorBeta$parameters[["components"]]
      betaDf <- 3.0
      sigmaDf <- prior$priorSigma$parameters[["df"]]
      sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			8
		},

		"tmix:non-informative" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["props"]]
      beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			9
		},

		"tmix:invChisq" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
      betaCov <- prior$priorBeta$parameters[["sigma"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["props"]]
      beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			10
		},

		"t:non-informative" = {
			betaMean <- prior$priorBeta@parameters[["mu"]]
			betaCov <- prior$priorBeta@parameters[["sigma"]]
			betaDf <- prior$priorBeta@parameters[["df"]]
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			2
		},

		"t:invChisq" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			2
		}
	)

	if(is.function(betaMean))
	  betaMean <- betaMean(dim.Cov)

	if(is.name(betaMean))
	  betaMean <- eval(call(betaMean, p = dim.Cov))

	if(is.function(betaCov))
	  betaCov <- betaCov( dim.Cov )

	if(is.name(betaCov))
    betaCov <- eval( call( betaCov, p = dim.Cov ) )

	if(is.vector(betaCov) && length(betaCov) == dim.Cov && dim.Cov > 1)
	  betaCov <- diag( betaCov )

  # the "normal:invChisq" && prior$conjugate  case
  if(iType == 6) 
	  betaCov <- betaCov / prior$priorBeta$parameters[["k0"]]

	# the mixture case
  if(is.element(iType, c( 7, 8, 9, 10))) {

    ##the t mixture case
    if(is.element(iType, c(9, 10))) {

      if(length(betaDf) == 1) {
        if(beta.components > 1)
          betaDf <- rep( betaDf, beta.components )

        else if(length(betaDf) != beta.components)
          stop("blm: degrees of freedom for all t mixture components must ",
               "be provided")
      }

      if(beta.components == length(beta.props) + 1)
        beta.props <- c(beta.props, 1 - sum(beta.props))

      if(any(beta.props < 0))
        stop("blm: problem with mixture components proportions: mixture ",
             "prior for the coefficients is not valid" )

      if(beta.k > 0 && beta.components == 2) {

        #number of mixtures is 2
        betaCov2 <- betaCov * (beta.k^2)
        betaCov <- cbind( betaCov, betaCov2 )

        if(is.vector(betaMean)) {

          #only one mean supplied
          if(length( betaMean) == dim.Cov)
                betaMean <- cbind( betaMean, betaMean )

          else if(dim.Cov != 1 || length(betaMean) != 2)
            stop("blm: problem with prior mean: mixture prior for the ",
                 "coefficients is not valid")
        }

        else if(is.matrix(betaMean)) {
          if(ncol(betaMean) != 2 || nrow(betaMean) != dim.Cov)
            stop("blm: problem with beta prior mean: mixture prior for ", 
                 "the coefficients is not valid")
        }
      }

      else if(beta.components > 2) {

        if(is.list(betaMean) && length(betaMean) == beta.components) {
          this.betaMean <- betaMean[[1]]

          for(i in seq(2, beta.components))
            this.betaMean <- cbind(this.betaMean, betaMean[[i]])

          betaMean <- this.betaMean
        }

        else if(is.vector(betaMean))
          betaMean <- rep(betaMean, beta.components)

        else if(is.matrix(betaMean)) {
          if(ncol(betaMean) != beta.components || nrow(betaMean) != dim.Cov)
            stop("blm: problem with beta prior mean: mixture prior for ",
                 "the coefficients is not valid")
        }

        else
          stop("blm: problem with beta prior mean: mixture prior for the ",
               "coefficients is not valid.")
      }

      if(is.list(betaCov) && length(betaCov) == beta.components) {

        this.betaCov <- NULL

        for(i in seq(1, beta.components)) {
          thisCov <- betaCov[[i]]

          if(is.vector(thisCov) && length(thisCov) == dim.Cov && dim.Cov > 1)
            thisCov <- diag(thisCov)

          else if(!is.matrix(thisCov))
            stop("blm: problem with prior covariance: mixture prior for ",
                 "the coefficients is not valid")

          if(i == 1)
            this.betaCov <- thisCov

          else
            this.betaCov <- cbind(this.betaCov, thisCov)
        }#end for i

        betaCov <- this.betaCov
      }

      else if(is.matrix(betaCov))
        betaCov <- rep( betaCov, beta.components )

      else if(is.vector(betaCov) && length(betaCov) == dim.Cov && dim.Cov > 1) {
        betaCov <- diag(betaCov)
        betaCov <- rep(betaCov, beta.components)    
      }

      else
        stop("blm: problem with prior covariance: mixture prior for ",
             "the coefficients is not valid")
    }
  }#end if iType is mixture

	## Call the C, C++ code ##
	generateInitPoints.blm(number.draws, X, Y, sigmaDf, sigmaScale, betaMean,
                         betaCov, betaDf, beta.components, beta.props)
}


# return: draws from beta and sigma2 based on mixture of prior and MLE's
#         (useful for creating different chains starting at the draws)

generateInitPoints.blm <- function(number.draws, X, Y, sigmaDf, sigmaScale,
                                   betaMean, betaCov, betaDf, beta.components,
                                   beta.props, mix.with.MLE = 1)
{
  ## mix.with.MLE
  # 1 => mix prior + MLE
  # 0 => prior
  # 2 => MLE

  number.data <- dim(X)[1]
  dim.Cov <- dim(X)[2]

  output.points <- matrix(0.0, nrow = dim.Cov + 1, ncol = number.draws)
  storage.mode(output.points) <- "double"

  mixture.MLE <- mix.with.MLE

  points <- .C("getInitialPointsBayesLM",
               as.integer(number.data),
               as.integer(dim.Cov),
               as.double(Y),
               as.double(t(X)),
               as.double(sigmaDf),
               as.double(sigmaScale),
               as.double(betaMean),
               as.double(betaCov),
               as.double(betaDf),
               as.double(beta.props),
               as.integer(beta.components),
               as.integer(number.draws),
               as.integer(mixture.MLE),
               output.points = output.points)

  #output is organized so that the first elements of it correspond to the samples
  #from beta (each beta at the time).
  #last number.draws elements contain draws of sigma2

  beta.draws <- matrix(points$output.points[seq(1, dim.Cov*number.draws)],
                       ncol = number.draws)

  sigma2.draws <- as.vector(points$output.points[seq(
                                                 dim.Cov * number.draws + 1,
                                                 (dim.Cov + 1)* number.draws)])

  list(beta = beta.draws, sigma = sigma2.draws)
}


