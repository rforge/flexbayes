blm <- function(formula, data, subset, weights, na.action, prior = blm.prior(),
                likelihood = blm.likelihood(), sampler = blm.sampler(),
                contrasts = NULL)
{
  ### if we want the statistics on variable drawings from Gibbs 
  ### or a mixture set print.stats to 1; otherwise set it to 0
  print.stats <- 1

  if(!missing(weights))
    stop("weighted regression is not implemented yet")

  #make sure formula is of the form y ~ 1  
  #  if(length(terms(formula)@term.labels)== 0 && terms(formula)@intercept == 1 &&  terms(formula)@response == 1)
  #  {
  #    location.scale <- T
  #  }
  #  else
  #  {
  #    location.scale <- F
  #  }

  cl <- mf <- match.call()
  keep <- c(TRUE, names(mf)[-1] %in% names(formals(model.frame.default)))
  mf <- mf[keep]
  mf[[1L]] <- quote(stats::model.frame)
  

#  terms.attrs <- attributes(terms(formula))
#  location.scale <- (!length(terms.attrs$term.labels) && terms.attrs$intercept && terms.attrs$response)

  #Kjell: these checks not necessary
  ###check if this is a call for the Location/Scale model
  #if(location.scale)
  #{
    #make sure formula is of the form y ~ 1  
    #if(length(terms(formula)@term.labels)!= 0)
    #{
    #  stop("blm: formula must contain no variables besides the response one.\n")
    #}

    #if(terms(formula)@intercept != 1)
    #{
    #  stop("blm: formula must contain the intercept term.")
    #}

    #if(terms(formula)@response == 0)
    #{
    #  stop("blm: must specified response variable in formula.")
    #}
  #}#end if location.scale


#model.blm <- call("model.frame", formula = formula, na.action = na.action)
#  if(!missing(data))
#   model.blm$data <- data

  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.extract(mf, "response")
  X <- model.matrix(mt, mf, contrasts)

  a <- attributes(mt)
  location.scale <- (!length(a$term.labels) && a$intercept && a$response)

  #Kjell: redundant!?
  #contrasts <- contrasts(X)

  n <- nrow(X)
  p <- ncol(X)

  ### First validate/interpret the arguments for the likelihood ###

  likelihood$type <- match.arg(likelihood$type, choices = c("normal", "t"))

  if(prior$conjugate && likelihood$type == "t")
    stop("the conjugate prior requires a normal likelihood")

  if(likelihood$df <= 0.0)
    stop("the degrees of freedom for the multivariate-t likelihood must be ",
         "positive")

  if(is.null(likelihood$errorCov))
    likelihood$errorCov <- 1

  else if(is.vector(likelihood$errorCov)) {
    if(length(likelihood$errorCov) == 1) {
      if(likelihood$errorCov <= 0) { 
        stop("the error covariance matrix must be positive")
      }
    }

    else {
      if(any(likelihood$errorCov)<= 0)
        stop ("the error covariance matrix must be positive.")

      if(length(likelihood$errorCov)!= n)
        stop("the dimension of the error covariance matrix must equal the number of observations")
    }
  }

  else if(is.matrix(likelihood$errorCov)) 
  {
    if(any(dim(likelihood$errorCov) != n))
      stop("the error covariance matrix must be square with dimension equal",
           "to the number of observations.")

    if (any(likelihood$errorCov - t(likelihood$errorCov) != 0))
      stop("the error covariance matrix is not symmetric")
  }

  ### Validate/interpret the arguments specifying the prior ###

  beta.props <- 1
  beta.components <- 1

#  if(class(prior$priorBeta) == "fbprior")
#    prior.beta <- prior$priorBeta$name
#  else
#    prior.beta <- "nonInformative"
#
#  if(class(prior$priorSigma) == "fbprior")
#    prior.sigma <- prior$priorSigma$name
#  else
#    prior.sigma <- "nonInformative"

  iType <- paste(likelihood$type,
                 prior$priorBeta$name,
                 prior$priorSigma$name,
                 sep = ":")

  iType <- switch(iType,

		"normal:nonInformative:nonInformative" = {
			betaMean <- rep(0.0, p)
			betaCov <- diag(p)
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			0
		},

		"normal:nonInformative:invChisq" = {
			betaMean <- rep(0.0, p)
			betaCov <- diag(p)
			betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			0
		},

		"normal:normal:nonInformative" = {
			betaMean <- prior$priorBeta$parameters$mean
			betaCov <- prior$priorBeta$parameters$S
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			1
		},

		"normal:normal:invChisq" = {
      betaMean <- prior$priorBeta$parameters$mean
      betaCov <- prior$priorBeta$parameters$S
      betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			ifelse(prior$conjugate, 6, 1)
		},

		"normal:t:nonInformative" = {
			betaMean <- prior$priorBeta$parameters[["mean"]]
			betaCov <- prior$priorBeta$parameters[["S"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			2
		},

		"normal:t:invChisq" = {
			betaMean <- prior$priorBeta$parameters[["mean"]]
			betaCov <- prior$priorBeta$parameters[["S"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			2
		},

		"normal:normalMixture:nonInformative" = {
			betaMean <- prior$priorBeta$parameters$mean
			betaCov <- prior$priorBeta$parameters$S
      beta.props <- prior$priorBeta$parameters$w
      beta.k <- prior$priorBeta$parameters$k0
      beta.components <- length(beta.props)
      #beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			7
		},

		"normal:normalMixture:invChisq" = {
      betaMean <- prior$priorBeta$parameters$mean
      betaCov <- prior$priorBeta$parameters$S
      beta.props <- prior$priorBeta$parameters$w
      beta.k <- prior$priorBeta$parameters$k0
      beta.components <- length(beta.props)
      #beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			ifelse(prior$conjugate, 15, 8)
		},

		"normal:tMixture:nonInformative" = {
			betaMean <- prior$priorBeta$parameters[["mean"]]
			betaCov <- prior$priorBeta$parameters[["S"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["w"]]
      beta.components <- length(beta.props)
      #beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			9
		},

		"normal:tMixture:invChisq" = {
			betaMean <- prior$priorBeta$parameters[["mean"]]
      betaCov <- prior$priorBeta$parameters[["S"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["w"]]
      beta.components <- length(beta.props)
      #beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			10
		},

		"t:nonInformative:nonInformative" = {
			betaMean <- rep(0.0, p)
			betaCov <- diag(p)
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			3
		},

		"t:nonInformative:invChisq" = {
			betaMean <- rep(0.0, p)
			betaCov <- diag(p)
			betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			3
		},

		"t:normal:nonInformative" = {
      betaMean <- prior$priorBeta$parameters$mean
      betaCov <- prior$priorBeta$parameters$S
      betaDf <- 3.0
      sigmaDf <- 0.0
      sigmaScale <- 1.0
			4
		},

		"t:normal:invChisq" = {
      betaMean <- prior$priorBeta$parameters$mean
      betaCov <- prior$priorBeta$parameters$S
      betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			4
		},

		"t:t:nonInformative" = {
			betaMean <- prior$priorBeta$parameters[["mean"]]
			betaCov <- prior$priorBeta$parameters[["S"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			5
		},

		"t:t:invChisq" = {
      betaMean <- prior$priorBeta$parameters[["mean"]]
      betaCov <- prior$priorBeta$parameters[["S"]]
      betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			5
		},

		"t:normalMixture:nonInformative" = {
      betaMean <- prior$priorBeta$parameters$mean
      betaCov <- prior$priorBeta$parameters$S
      beta.props <- prior$priorBeta$parameters$w
      beta.k <- prior$priorBeta$parameters$k0
      beta.components <- length(beta.props)
      #beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			11
		},

		"t:normalMixture:invChisq" = {
      betaMean <- prior$priorBeta$parameters$mean
      betaCov <- prior$priorBeta$parameters$S
      beta.props <- prior$priorBeta$parameters$w
      beta.k <- prior$priorBeta$parameters$k0
      beta.components <- length(beta.props)
      #beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			12
		},

		"t:tMixture:nonInformative" = {
      betaMean <- prior$priorBeta$parameters[["mean"]]
      betaCov <- prior$priorBeta$parameters[["S"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["w"]]
      beta.components <- length(beta.props)
      #beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			13
		},

		"t:tMixture:invChisq" = {
      betaMean <- prior$priorBeta$parameters[["mean"]]
      betaCov <- prior$priorBeta$parameters[["S"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["w"]]
      beta.components <- length(beta.props)
      #beta.components <- prior$priorBeta$parameters[["components"]]
      betaDf <- prior$priorBeta$parameters[["df"]]
      sigmaDf <- prior$priorSigma$parameters[["df"]]
      sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			14
		}
	)


#  if(is.function(betaMean))
#    betaMean <- betaMean(p)

#  if(is.name(betaMean))
#    betaMean <- eval(call(betaMean, p))

#  if(is.function(betaCov))
#    betaCov <- betaCov(p)

#  if(is.name(betaCov))
#    betaCov <- eval(call(betaCov, p))

#  if(is.vector(betaCov)) {
#    if(length(betaCov) == p && p > 1)
#      betaCov <- diag(betaCov)
#    else if(length(betaCov) == 1 && p > 1)
#      betaCov <- betaCov * diag(p)
#  }

  # the mixture case
  if(iType %in% c(7, 8, 9, 10, 11, 12, 13, 14, 15)) {
    ##the t mixture case
    if(is.element(iType, c(9, 10, 13, 14)))
    {
      if(length(betaDf) == 1) {
        if(beta.components > 1)
          betaDf <- rep(betaDf, beta.components)
      }
      else if(length(betaDf) != beta.components)
        stop("blm: degrees of freedom for all t mixture components must be provided.\n")
    }#end t-mixture case

    if(beta.components == length(beta.props) + 1)
      beta.props <- c(beta.props, 1 - sum(beta.props))

    else if(length(beta.props)< beta.components) {

      if(any(beta.props < 0))
        stop("blm: problem with mixture components proportions. mixture prior for the coefficients is not valid.")

      if(sum(beta.props) >= 1)
        stop("blm: proportions for all mixture components must be provided.\n")

      else {
        prop.left <- 1.0 - sum(beta.props)
        props.left <- rep(prop.left /(beta.components - length(beta.props)), beta.components - length(beta.props))
        beta.props <- c(beta.props, props.left)
        warning("blm: proportions for all mixture components have been automatically provided.\n")
      }
    }

    if(any(beta.props < 0))
      stop("blm: problem with mixture components proportions. mixture prior for the coefficients is not valid.")

    if(beta.components > 1) {
      if(is.list(betaMean) && length(betaMean) == beta.components) {
        this.betaMean <- betaMean[[1]]

        for(i in seq(2, beta.components))
          this.betaMean <- cbind(this.betaMean, betaMean[[ i ]], deparse.level = 0)

        betaMean <- this.betaMean
      }

      else if(is.vector(betaMean)) {

        if(length(betaMean) == p)
          betaMean <- matrix(rep(betaMean, beta.components), ncol = beta.components)

        else if(length(betaMean)== beta.components && p == 1)
          betaMean <- matrix(betaMean, ncol = beta.components)

        else
          stop("blm: problem with beta prior mean. mixture prior for the coefficients is not valid.")
      }

      else if(is.matrix(betaMean)) {
        if(ncol(betaMean) != beta.components || nrow(betaMean) != p)
          stop("blm: problem with beta prior mean. mixture prior for the coefficients is not valid.")
      }

      else
        stop("blm: problem with beta prior mean. mixture prior for the coefficients is not valid.")

      if(is.list(betaCov) && length(betaCov) == beta.components) {
        this.betaCov <- NULL
        for(i in seq(1:beta.components)) {
          thisCov <- betaCov[[i]]

          if(is.vector(thisCov) && length(thisCov) == p && p > 1)
            thisCov <- diag(thisCov)

          else if(is.vector(thisCov) && length(thisCov) == 1) {
            if(p > 1)
              thisCov <- diag(rep(thisCov, p))
          }

          else if(!is.matrix(thisCov))
            stop("blm: problem with prior covariance. mixture prior for the coefficients is not valid.")

          if(i == 1)
            this.betaCov <- thisCov

          else
            this.betaCov <- cbind(this.betaCov, thisCov)
        }#end for i

        betaCov <- this.betaCov
      }

      else if(is.matrix(betaCov))
        betaCov <- matrix(rep(betaCov, beta.components), ncol =  p * beta.components)

      else if(is.vector(betaCov) && length(betaCov) == p) {

        if(p > 1)
          betaCov <- diag(betaCov)

        betaCov <- matrix(rep(betaCov, beta.components), ncol =  p * beta.components)  
      }

      else {
        stop("blm: problem with prior covariance. mixture prior for the coefficients is not valid.")
      }
    }

    else
      stop("blm: number of components in a prior mixture must be larger than one.")
  }

  ###  THIS NEEDS MORE REVISION        ###
  ###  CHANGE GUI TO ASK FOR VARIANCE  ###

  ## For the Location Scale Model betaCov corresponds to standard deviations
  if(location.scale && p == 1)
      betaCov <- betaCov^2

  # the "normal:normal:invChisq" && prior$conjugate  case
  if(iType == 6 || iType == 15)
    betaCov <- betaCov / prior$priorBeta$parameters[["k0"]]

  if(prior$conjugate && !is.element(iType, c(6, 15))) {
    prior$conjugate <- FALSE
    warning("blm: prior has been set up as conjugate, but provided priors",
            "do not match the conjugate type; prior has been modified to",
            "non-conjugate type")
  }

  # the mixture case
  if(is.element(iType, c(7, 8, 9, 10, 11, 12, 13, 14, 15)))
  {
    ##now look at covariance/variance inflation factors
    inflation.factors <- beta.k
    if(length(beta.k) > 0 && all(beta.k > 0) && any(beta.k != 1) && beta.components >= 2) {

      if(length(beta.k) < beta.components - 1)
        stop("blm: too few covariance inflation factor(s) provided.")

      else if(length(beta.k) > beta.components)
        stop("blm: too many covariance inflation factor(s) provided.")

      ##make length of beta.k be beta.components
      if(length(beta.k) == beta.components - 1)
        beta.k <- c(1, beta.k)

      inflation.factors <- rep(beta.k[1], p)
      for(i in seq(2, beta.components))
        inflation.factors <- c(inflation.factors, rep(beta.k[i], p))
    
      inflation.factors <- diag(inflation.factors^2)
      betaCov <- betaCov %*% inflation.factors
    }

    if(beta.components >= 2) {
      ### check if all components in mixture are equal

      idx.s <- 1
      idx.e <- p

      for(i in seq(1, beta.components - 1)) {
        first.mean <- betaMean[, i]
        first.Cov <- betaCov[, seq(idx.s, idx.e)]
        jdx.s <- idx.e + 1
        jdx.e <- idx.e + p
        for(j in seq(i + 1, beta.components)) {
          second.mean <- betaMean[, j]
          second.Cov <- betaCov[, seq(jdx.s, jdx.e)]
          diff <- max(abs(first.mean - second.mean)) + max(abs(first.Cov - second.Cov))

          if(diff == 0)
            stop("blm: Mixture components ", i, " and ", j, " coincide.")

          jdx.s <- jdx.e + 1
          jdx.e <- jdx.e + p
        }#end for j

        idx.s <- idx.e + 1
        idx.e <- idx.e + p
      }#end for i
    }

  }#end if mixture

  # the non-mixture case
  if(!is.element(iType, c(7, 8, 9, 10, 11, 12, 13, 14, 15))) {

    if(length(betaMean) != p) {
      stop("blm: wrong dimension for coefficients mean prior; it is ",
            length(betaMean), " and should be ", p)
    }
  }

  
  #Kjell: lets not do this anymore
  ### set the seed

  #set.seed(random.seed)

  ### Get the control arguments ###

  #burnInLength <- sampler$control$bSize
  #simulationsToPerform <- sampler$control$simSize
  #sampleFrequency <- sampler$control$freqSize
  
  burnInLength <- sampler$control$nBurnin
  simulationsToPerform <- sampler$control$nSamples
  sampleFrequency <- sampler$control$nThin
  
  ## Get the starting point for the simulation ###

  nChains <- sampler$nChains
  init.point <- sampler$start

  if(is.null(sampler$sampler))
    sampler.type <- 0

  else if(sampler$sampler == "Gibbs")
    sampler.type <- 0

  else #exact
    sampler.type <- 1

  if(is.null(init.point) || init.point$type == "prior") {

    # flag for prior means starting points

    ## flag for specify/mixture starting points
    read.init.point <- 1

      #warning("blm: need ", nChains, " initial starting points for simulations.\n",
      #"Drawing points at random from prior.")

      # draw nChains initial values at random from prior

    starting.points <- generateInitPoints.blm(nChains, X, Y, sigmaDf, sigmaScale, 
                                              betaMean, betaCov, betaDf, beta.components,
                                              beta.props, mix.with.MLE = 0)
  }

  else {
    # flag for specify/mixture starting points
    read.init.point <- 1

    if(init.point$type == "user's choice") {
      if(is.vector(init.point$beta)) {
        if(length(init.point$beta) == p) {
          if(nChains == 1)
            beta.init <- init.point$beta

          else
            stop("blm: need to specify ", nChains, " starting vectors for beta")
        }

        else if(length(init.point$beta) != nChains && p == 1)
          stop("blm: need to specify ", nChains, " starting values for beta")

        else if(p == 1)
          beta.init <- init.point$beta

        else
          stop("dimension of initial beta [", length(init.point$beta),
                "] does not match dimension of parameter beta [", p, "]")

        beta.init <- matrix(beta.init, ncol = nChains)

      }#end if beta is vector

      else if(is.matrix(init.point$beta)) {

        if(ncol(init.point$beta) != nChains)
          stop("blm: need to specify ", nChains, " starting vectors for beta.")

        else if(nrow(init.point$beta) != p)
          stop("dimension of initial beta [", nrow(init.point$beta),
                "] does not match dimension of parameter beta [", p, "]")

        else
          beta.init <- init.point$beta
      }

      if(is.vector(init.point$sigma)) {
        if(length(init.point$sigma) != nChains)
          stop("blm: need to specify ", nChains, " starting values for scale")

        else if(any(init.point$sigma <= 0))
          stop("blm: starting values for scale must be positive")

        else
          sigma.init <- (init.point$sigma)^2
      }

      else
        stop("blm: wrong starting values argument for scale")

      starting.points <- list(beta = beta.init, sigma = sigma.init)

    }#end if "specify"

    else if(init.point$type == "prior + likelihood") {
      # draw nChains initial values at random from MLE-prior mixture
      starting.points <- generateInitPoints.blm(nChains, X, Y, sigmaDf, sigmaScale, 
                                                betaMean, betaCov, betaDf, beta.components,
                                                beta.props)
    }

    else if(init.point$type == "likelihood") {
      # draw nChains initial values at random from MLE
      starting.points <- generateInitPoints.blm(nChains, X, Y, sigmaDf, sigmaScale, 
                                                betaMean, betaCov, betaDf, beta.components,
                                                beta.props, mix.with.MLE = 2)
    }

  }#end if not null(init.point)


  ### Call C++ code ###


  ##update prior parameters
#  if(class(prior$priorBeta) == "fbprior")
#  {
#    prior$priorBeta$parameters[["mu"]] <- betaMean
#    prior$priorBeta$parameters[["sigma"]] <- betaCov
#
#    if(is.element(prior$priorBeta$name, c("normmix", "tmix"))) {
#      prior$priorBeta$parameters[["props"]] <- beta.props
#
#      if(prior$priorBeta$name == "tmix")
#        prior$priorBeta$parameters[["df"]] <- betaDf
#    }
#  }

#  else #create proper prior
#    prior$priorBeta <- bayes.nonInformative(betaMean, betaCov)
#
#  if(class(prior$priorSigma)== "fbprior") {
#    prior$priorSigma$parameters[["df"]] <- sigmaDf
#    prior$priorSigma$parameters[["sigma0.sq"]] <- sigmaScale
#  }

#  else #create proper prior
#    prior$priorSigma <- bayes.invChisq(sigmaDf, sigmaScale)

#if(location.scale && length(betaMean) == 1) #change name of "intercept" to "mean"
#    beta.names <- "MEAN"

  #run the code
  coefs <- list()
  sigma <- list()
  tau2ErrorScale <- list()

  for(i in 1:nChains) {

    fit <- blm.fit(X,
                   Y,
                   likelihood$errorCov,
                   likelihood$df,
                   sigmaDf,
                   sigmaScale,
                   iType,
                   betaMean,
                   betaCov,
                   betaDf,
                   beta.props,
                   beta.components,
                   burnInLength,
                   simulationsToPerform,
                   sampleFrequency,
                   sampler.type,
                   read.init.point,
                   starting.points$beta[, i],
                   starting.points$sigma[i],
                   print.stats)

    #create new object blm
#    gibbs.drawing.stats <- NULL
#    mixture.drawing.stats <- NULL

#    if(length(fit$gibbs) > 1) {
#      sum.gibbs <- sum(fit$gibbs)
#
#      if(sum.gibbs > 0) {
#        gibbs.drawing.stats <- matrix(c(fit$gibbs / sum.gibbs, fit$gibbs), ncol = 2)
#        gibbs.names <- c("beta", "sigma")
#
#        if(iType %in% c(2, 5))
#          gibbs.names <- c(gibbs.names, "TAU:BETA")
#
#        if(iType %in% c(3, 4, 5, 11, 12, 13, 14))
#          gibbs.names <- c(gibbs.names, paste("TAU:ERROR:", seq(1, n), sep = ""))
#
#        dimnames(gibbs.drawing.stats) <- list(gibbs.names, c("Proportion", "Count"))
#      }
#    }

#    if(length(fit$mixture) > 1) {
#      sum.mixture <- sum(fit$mixture)
#
#      if(sum.mixture > 0) {
#        mixture.drawing.stats <- matrix(c(fit$mixture / sum.mixture, fit$mixture), ncol = 2)
#        dimnames(mixture.drawing.stats)<- list(paste("Component:", seq(1, beta.components), sep = ""),
#                                                   c("Proportion", "Count"))
#      }
#    }

#	  if(iType %in% c(3, 4, 5, 11, 12, 13, 14)) {
#    	tau2ErrorScale.names <- paste("TAU:ERROR:", seq(1, dim(X)[1]), sep = "")
#    	samples.this.chain <- fit$samples[, 1:(p + 1 + dim(X)[1])]
#    	dimnames(samples.this.chain) <- list(NULL, c(beta.names, "~sigma~", tau2ErrorScale.names))
#
#  		blmodel[[i]] <- mcmc(samples.this.chain,
#                           thin = sampleFrequency,
#  		                     start = burnInLength + 1,
#                           end = simulationsToPerform*sampleFrequency + burnInLength)
#    }
#
#    else {
#    	samples.this.chain <- fit$samples[, 1:(p + 1)]
#    	dimnames(samples.this.chain) <- list(NULL, c(beta.names, "~sigma~"))
#
#  		blmodel[[i]] <- mcmc(samples.this.chain,
#                           thin = sampleFrequency,
#                           start = burnInLength + 1,
#                           end = simulationsToPerform*sampleFrequency + burnInLength)
#    }

    samples.this.chain <- fit$samples[, 1:(p + 1)]
    dimnames(samples.this.chain) <- list(NULL, c(dimnames(X)[[2]], "~sigma~"))

    samples <- fit$samples[, 1:p, drop = FALSE]
    dimnames(samples) <- list(NULL, colnames(X))

    coefs[[i]] <- mcmc(samples,
                       thin = sampleFrequency,
                       start = burnInLength + 1,
                       end = simulationsToPerform*sampleFrequency + burnInLength)

    samples <- fit$samples[, p + 1, drop = FALSE]
    dimnames(samples) <- list(NULL, "sigma")

    sigma[[i]] <- mcmc(samples,
                       thin = sampleFrequency,
                       start = burnInLength + 1,
                       end = simulationsToPerform*sampleFrequency + burnInLength)

    if((d <- dim(fit$samples)[2]) > p + 1) {

      samples <- fit$samples[, (p + 2):d, drop = FALSE]
      dimnames(samples) <- list(NULL, paste("TAU:ERROR:", seq(1, n), sep = ""))

      tau2ErrorScale[[i]] <- mcmc(samples,
                                  thin = sampleFrequency,
                                  start = burnInLength + 1,
                                  end = simulationsToPerform*sampleFrequency + burnInLength)

    }
  }

	ans <- list(coefficients = mcmc.list(coefs),
              sigma = mcmc.list(sigma),
              tau2ErrorScale = mcmc.list(tau2ErrorScale),
              call = cl,
              terms = mt,
              model = mf,
              contrasts = contrasts)

  oldClass(ans) <- "blm"
  ans
}



