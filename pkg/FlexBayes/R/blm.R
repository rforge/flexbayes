blm <- function(formula, data, prior = blm.prior(), likelihood = blm.likelihood(),
                sampler = blm.sampler(), random.seed = .Random.seed, na.action = na.fail,
                contrasts = NULL)
{
  ### if we want the statistics on variable drawings from Gibbs 
  ### or a mixture set print.stats to 1; otherwise set it to 0
  print.stats <- 1

  #make sure formula is of the form y ~ 1  
  #  if(length(terms(formula)@term.labels)== 0 && terms(formula)@intercept == 1 &&  terms(formula)@response == 1)
  #  {
  #    location.scale <- T
  #  }
  #  else
  #  {
  #    location.scale <- F
  #  }

  terms.attrs <- attributes(terms(formula))
  location.scale <- (!length(terms.attrs$term.labels) && terms.attrs$intercept && terms.attrs$response)

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

  model.blm <- call("model.frame", formula = formula, na.action = na.action)
  if(!missing(data))
    model.blm$data <- data

  model.blm <- eval(model.blm, parent.frame())
  Terms <- attr(model.blm, "terms")
  Y <- model.extract(model.blm, "response")
  X <- model.matrix(Terms, model.blm, contrasts)

  #Kjell: redundant!?
  #contrasts <- contrasts(X)

  dim.Cov <- dim(X)[2]
  number.data <- length(Y)

  ### First validate/interpret the arguments for the likelihood ###

  degreesOfFreedom.likelihood <- likelihood$df
  likelihood$type <- match.arg(likelihood$type, choices = c("normal", "t"))

  if(prior$conjugate && likelihood$type == "t")
    stop("The conjugate prior requires a normal likelihood.")

  if(degreesOfFreedom.likelihood <= 0.0) {
    warning(paste("the degrees of freedom for the multivariate-t likelihood must be",
	                "positive; using the default value of 3", sep = "\n"))
    degreesOfFreedom.likelihood <- 3
  }

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

      if(length(likelihood$errorCov)!= number.data)
        stop("the dimension of the error covariance matrix must equal the number of observations")
    }
  }

  else if(is.matrix(likelihood$errorCov)) 
  {
    if(any(dim(likelihood$errorCov) != number.data))
      stop("the error covariance matrix must be square with dimension equal",
           "to the number of observations.")

    if (any(likelihood$errorCov - t(likelihood$errorCov) != 0))
      stop("the error covariance matrix is not symmetric")
  }

  ### Validate/interpret the arguments specifying the prior ###

  beta.props <- 1
  beta.components <- 1

  if(class(prior$priorBeta) == "fbdstn")
    prior.beta <- prior$priorBeta$name
  else
    prior.beta <- "non-informative"

  if(class(prior$priorSigma) == "fbdstn")
    prior.sigma <- prior$priorSigma$name
  else
    prior.sigma <- "non-informative"

  iType <- paste(likelihood$type, prior.beta, prior.sigma, sep = ":")

  iType <- switch(iType,

		"normal:non-informative:non-informative" = {
			betaMean <- rep(0.0, dim.Cov)
			betaCov <- diag(dim.Cov)
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			0
		},

		"normal:non-informative:invChisq" = {
			betaMean <- rep(0.0, dim.Cov)
			betaCov <- diag(dim.Cov)
			betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			0
		},

		"normal:normal:non-informative" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			1
		},

		"normal:normal:invChisq" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]

      if(prior$conjugate && dim.Cov == 1)
        betaCov <- diag(dim.Cov)
      else
        betaCov <- prior$priorBeta$parameters[["sigma"]]

      ##after evaluating betaCov (it could be a function), divide by k0 (see below) 
			betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			ifelse(prior$conjugate, 6, 1)
		},

		"normal:t:non-informative" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			2
		},

		"normal:t:invChisq" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			2
		},

		"normal:normal mixture:non-informative" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["props"]]
      beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			7
		},

		"normal:normal mixture:invChisq" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
      betaCov <- prior$priorBeta$parameters[["sigma"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["props"]]
      beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			ifelse(prior$conjugate, 15, 8)
		},

		"normal:t mixture:non-informative" = {
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

		"normal:t mixture:invChisq" = {
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

		"t:non-informative:non-informative" = {
			betaMean <- rep(0.0, dim.Cov)
			betaCov <- diag(dim.Cov)
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			3
		},

		"t:non-informative:invChisq" = {
			betaMean <- rep(0.0, dim.Cov)
			betaCov <- diag(dim.Cov)
			betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			3
		},

		"t:normal:non-informative" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			4
		},

		"t:normal:invChisq" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
			betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			4
		},

		"t:t:non-informative" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["igma"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			5
		},

		"t:t:invChisq" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			5
		},

		"t:normal mixture:non-informative" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["props"]]
      beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- 3.0
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			11
		},

		"t:normal mixture:invChisq" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
      betaCov <- prior$priorBeta$parameters[["sigma"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["props"]]
      beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- 3.0
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			12
		},

		"t:t mixture:non-informative" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
			betaCov <- prior$priorBeta$parameters[["sigma"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["props"]]
      beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- 0.0
			sigmaScale <- 1.0
			13
		},

		"t:t mixture:invChisq" = {
			betaMean <- prior$priorBeta$parameters[["mu"]]
      betaCov <- prior$priorBeta$parameters[["sigma"]]
      beta.k <- prior$priorBeta$parameters[["k"]]
      beta.props <- prior$priorBeta$parameters[["props"]]
      beta.components <- prior$priorBeta$parameters[["components"]]
			betaDf <- prior$priorBeta$parameters[["df"]]
			sigmaDf <- prior$priorSigma$parameters[["df"]]
			sigmaScale <- prior$priorSigma$parameters[["sigma0.sq"]]
			14
		}
	)


  if(is.function(betaMean))
    betaMean <- betaMean(dim.Cov)

  if(is.name(betaMean))
    betaMean <- eval(call(betaMean, p = dim.Cov))

  if(is.function(betaCov))
    betaCov <- betaCov(dim.Cov)

  if(is.name(betaCov))
    betaCov <- eval(call(betaCov, p = dim.Cov))

  if(is.vector(betaCov)) {
    if(length(betaCov) == dim.Cov && dim.Cov > 1)
      betaCov <- diag(betaCov)
    else if(length(betaCov) == 1 && dim.Cov > 1)
      betaCov <- betaCov * diag(dim.Cov)
  }

  # the mixture case
  if(is.element(iType, c(7, 8, 9, 10, 11, 12, 13, 14, 15)))
  {
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

        if(length(betaMean) == dim.Cov)
          betaMean <- matrix(rep(betaMean, beta.components), ncol = beta.components)

        else if(length(betaMean)== beta.components && dim.Cov == 1)
          betaMean <- matrix(betaMean, ncol = beta.components)

        else
          stop("blm: problem with beta prior mean. mixture prior for the coefficients is not valid.")
      }

      else if(is.matrix(betaMean)) {
        if(ncol(betaMean) != beta.components || nrow(betaMean) != dim.Cov)
          stop("blm: problem with beta prior mean. mixture prior for the coefficients is not valid.")
      }

      else
        stop("blm: problem with beta prior mean. mixture prior for the coefficients is not valid.")

      if(is.list(betaCov) && length(betaCov) == beta.components) {
        this.betaCov <- NULL
        for(i in seq(1: beta.components)) {
          thisCov <- betaCov[[i]]

          if(is.vector(thisCov) && length(thisCov) == dim.Cov && dim.Cov > 1)
            thisCov <- diag(thisCov)

          else if(is.vector(thisCov) && length(thisCov) == 1) {
            if(dim.Cov > 1)
              thisCov <- diag(rep(thisCov, dim.Cov))
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
        betaCov <- matrix(rep(betaCov, beta.components), ncol =  dim.Cov * beta.components)

      else if(is.vector(betaCov) && length(betaCov) == dim.Cov) {

        if(dim.Cov > 1)
          betaCov <- diag(betaCov)

        betaCov <- matrix(rep(betaCov, beta.components), ncol =  dim.Cov * beta.components)  
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
  if(location.scale && dim.Cov == 1)
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

      inflation.factors <- rep(beta.k[1], dim.Cov)
      for(i in seq(2, beta.components))
        inflation.factors <- c(inflation.factors, rep(beta.k[i], dim.Cov))
    
      inflation.factors <- diag(inflation.factors^2)
      betaCov <- betaCov %*% inflation.factors
    }

    if(beta.components >= 2) {
      ### check if all components in mixture are equal

      idx.s <- 1
      idx.e <- dim.Cov

      for(i in seq(1, beta.components - 1)) {
        first.mean <- betaMean[, i]
        first.Cov <- betaCov[, seq(idx.s, idx.e)]
        jdx.s <- idx.e + 1
        jdx.e <- idx.e + dim.Cov
        for(j in seq(i + 1, beta.components)) {
          second.mean <- betaMean[, j]
          second.Cov <- betaCov[, seq(jdx.s, jdx.e)]
          diff <- max(abs(first.mean - second.mean)) + max(abs(first.Cov - second.Cov))

          if(diff == 0)
            stop("blm: Mixture components ", i, " and ", j, " coincide.")

          jdx.s <- jdx.e + 1
          jdx.e <- jdx.e + dim.Cov
        }#end for j

        idx.s <- idx.e + 1
        idx.e <- idx.e + dim.Cov
      }#end for i
    }

  }#end if mixture

  # the non-mixture case
  if(!is.element(iType, c(7, 8, 9, 10, 11, 12, 13, 14, 15))) {

    if(length(betaMean)!= dim.Cov) {
      stop("blm: wrong dimension for coefficients mean prior; it is ",
            length(betaMean), " and should be ", dim.Cov)
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

  number.chains <- sampler$number.chains
  init.point <- sampler$init.point

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

      #warning("blm: need ", number.chains, " initial starting points for simulations.\n",
      #"Drawing points at random from prior.")

      # draw number.chains initial values at random from prior

    starting.points <- generateInitPoints.blm(number.chains, X, Y, sigmaDf, sigmaScale, 
                                              betaMean, betaCov, betaDf, beta.components,
                                              beta.props, mix.with.MLE = 0)


    #}
  }

  else {
    # flag for specify/mixture starting points
    read.init.point <- 1

    if(init.point$type == "user's choice") {
      if(is.vector(init.point$beta)) {
        if(length(init.point$beta) == dim.Cov) {
          if(number.chains == 1)
            beta.init <- init.point$beta

          else
            stop("blm: need to specify ", number.chains, " starting vectors for beta")
        }

        else if(length(init.point$beta) != number.chains && dim.Cov == 1)
          stop("blm: need to specify ", number.chains, " starting values for beta")

        else if(dim.Cov == 1)
          beta.init <- init.point$beta

        else
          stop("dimension of initial beta [", length(init.point$beta),
                "] does not match dimension of parameter beta [", dim.Cov, "]")

        beta.init <- matrix(beta.init, ncol = number.chains)

      }#end if beta is vector

      else if(is.matrix(init.point$beta)) {

        if(ncol(init.point$beta) != number.chains)
          stop("blm: need to specify ", number.chains, " starting vectors for beta.")

        else if(nrow(init.point$beta) != dim.Cov)
          stop("dimension of initial beta [", nrow(init.point$beta),
                "] does not match dimension of parameter beta [", dim.Cov, "]")

        else
          beta.init <- init.point$beta
      }

      if(is.vector(init.point$sigma)) {
        if(length(init.point$sigma) != number.chains)
          stop("blm: need to specify ", number.chains, " starting values for scale")

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
      # draw number.chains initial values at random from MLE-prior mixture
      starting.points <- generateInitPoints.blm(number.chains, X, Y, sigmaDf, sigmaScale, 
                                                betaMean, betaCov, betaDf, beta.components,
                                                beta.props)
    }

    else if(init.point$type == "likelihood") {
      # draw number.chains initial values at random from MLE
      starting.points <- generateInitPoints.blm(number.chains, X, Y, sigmaDf, sigmaScale, 
                                                betaMean, betaCov, betaDf, beta.components,
                                                beta.props, mix.with.MLE = 2)
    }

  }#end if not null(init.point)


  ### Call C++ code ###

  #get variable names for output
  beta.names <- dimnames(X)[[2]]




  #Kjell: redundant?



  ##update prior parameters
  if(class(prior$priorBeta) == "fbdstn")
  {
    prior$priorBeta$parameters[["mu"]] <- betaMean
    prior$priorBeta$parameters[["sigma"]] <- betaCov

    if(is.element(prior$priorBeta$name, c("normmix", "tmix"))) {
      prior$priorBeta$parameters[["props"]] <- beta.props

      if(prior$priorBeta$name == "tmix")
        prior$priorBeta@parameters[["df"]] <- betaDf
    }
  }

  else #create proper prior
    prior$priorBeta <- bayes.nonInformative(betaMean, betaCov)

  if(class(prior$priorSigma)== "fbdstn") {
    prior$priorSigma@parameters[["df"]] <- sigmaDf
    prior$priorSigma@parameters[["sigma0.sq"]] <- sigmaScale
  }

  else #create proper prior
    prior$priorSigma <- bayes.invChisq(sigmaDf, sigmaScale)

  if(location.scale && length(betaMean) == 1) #change name of "intercept" to "mean"
    beta.names <- "MEAN"

  #run the code
  blmodel <- vector("list", number.chains)

  for(i in seq(1, number.chains)) {

    fit <- fit.bayeslm(X, Y,
                       likelihood$errorCov,
                       degreesOfFreedom.likelihood,
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
    gibbs.drawing.stats <- NULL
    mixture.drawing.stats <- NULL

    if(length(fit$gibbs) > 1) {
      sum.gibbs <- sum(fit$gibbs)

      if(sum.gibbs > 0) {
        gibbs.drawing.stats <- matrix(c(fit$gibbs / sum.gibbs, fit$gibbs), ncol = 2)
        gibbs.names <- c("beta", "sigma")

        if(is.element(iType, c(2, 5)))
          gibbs.names <- c(gibbs.names, "TAU:BETA")

        if(is.element(iType, c(3, 4, 5, 11, 12, 13, 14)))
          gibbs.names <- c(gibbs.names, paste("TAU:ERROR:", seq(1, number.data), sep = ""))

        dimnames(gibbs.drawing.stats) <- list(gibbs.names, c("Proportion", "Count"))
      }
    }
   
    if(length(fit$mixture) > 1) {
      sum.mixture <- sum(fit$mixture)

      if(sum.mixture > 0) {
        mixture.drawing.stats <- matrix(c(fit$mixture / sum.mixture, fit$mixture), ncol = 2)
        dimnames(mixture.drawing.stats)<- list(paste("Component:", seq(1, beta.components), sep = ""),
                                                   c("Proportion", "Count"))
      }
    }

	  if(is.element(iType, c(3, 4, 5, 11, 12, 13, 14))) {
    	tau2ErrorScale.names = paste("TAU:ERROR:", seq(1, dim(X)[1]), sep = "")
    	samples.this.chain <- fit$samples[, 1:(dim.Cov + 1 + dim(X)[1])]
    	dimnames(samples.this.chain) <- list(NULL, c(beta.names, "SIGMA", tau2ErrorScale.names))

  		blmodel[[i]] <- mcmc(samples.this.chain,
                           thin = sampleFrequency, 
  		                     start = burnInLength + 1,
                           end = simulationsToPerform*sampleFrequency + burnInLength)
    }

    else {
    	samples.this.chain <- fit$samples[, 1:(dim.Cov + 1)]
    	dimnames(samples.this.chain) <- list(NULL, c(beta.names, "SIGMA"))
    	
  		blmodel[[i]] <- mcmc(samples.this.chain,
                           thin = sampleFrequency, 
                           burnin = burnInLength)
    }
  }#end for loop

	posterior(sims = mcmc.list(blmodel), call = match.call())
}#end blm generator


fit.bayeslm <- function(X, Y, errorCov, degreesOfFreedom.likelihood, sigmaDf,
                        sigmaScale, iType, betaMean, betaCov, betaDf, props = 1,
                        n.mixtures = 1, burnInLength, simulationsToPerform,
                        sampleFrequency, sampler.type = 0, read.init.point = 0,
                        beta.init, sigma.init, print.stats = 0)
{

  if(n.mixtures == 1)
    props <- 1

  number.data <- dim(X)[1]
	dim.Cov <- dim(X)[2]
	dim.error.Cov <- length(errorCov)

  #cat("sampler type is ", sampler.type , "  iType is ", iType, "\n")

  keep_tau2_errors <- 0
  #account for beta and sigma
  number.vars <- 2 

  if(is.element(iType, c(2, 5))) #account for tau2 for beta t-prior
    number.vars <- number.vars + 1

  if(is.element(iType, c(3, 4, 5, 11, 12, 13, 14))) {

    ##t-likelihood case (keep tau2 errors in data augmentation)
	  output.samples <- matrix(0.0, nrow = simulationsToPerform, ncol = dim.Cov + 1 + number.data)
    keep_tau2_errors <- 1

    #account for all tau2_errors
    number.vars <- number.vars + number.data
  }

  else ##exact sampling
    output.samples <- matrix(0.0, nrow = simulationsToPerform, ncol = dim.Cov + 1)

  if(print.stats == 1) {

    #the Gibbs sampler case
    if(sampler.type == 0)
      gibbs.stats <- rep(0, number.vars)

    else
      gibbs.stats <- 0

    #the mixture prior case
    if(is.element(iType, c(7, 8, 9, 10, 11, 12, 13, 14, 15)))
      mixture.stats <- rep(0, n.mixtures)

    else
      mixture.stats <- 0
  }

	fit <- .C("fitBayesianLM",
            as.integer(number.data),
            as.integer(dim.Cov),
	          as.double(Y),
	          as.double(t(X)),
	          as.integer(dim.error.Cov),
	          as.double(errorCov),
	          as.double(degreesOfFreedom.likelihood),
	          as.double(sigmaDf),
	          as.double(sigmaScale),
	          as.integer(iType), 
	          as.double(betaMean),
	          as.double(betaCov),
	          as.double(betaDf),
	          as.double(props),
	          as.integer(n.mixtures),
	          as.integer(burnInLength), 
	          as.integer(simulationsToPerform),
	          as.integer(sampleFrequency),
	          as.integer(read.init.point),
	          as.double(beta.init),
	          as.double(sigma.init),
	          as.integer(sampler.type),
	          as.integer(print.stats),
	          output.samples = as.double(output.samples),
	          gibbs.stats = as.double(gibbs.stats),
	          mixture.stats = as.double(mixture.stats))

	# output is organized as a matrix of simulationsToPerform x beta dimension
	# the column beyond beta dimension is associated to simulations of sigma
	# (= sqrt(sigma2))
  # if keep tau2 errors then the last number.data columns contain the simulated
  # sqrt(scale)of each tau2_i

  if(keep_tau2_errors == 0)
    out.samples <- matrix(fit$output.samples, nrow = simulationsToPerform, ncol = dim.Cov + 1)

  else
    out.samples <- matrix(fit$output.samples, nrow = simulationsToPerform, ncol = dim.Cov + 1 + number.data)

  list(samples = out.samples,
       gibbs = fit$gibbs.stats,
       mixture = fit$mixture.stats)
}




################################################################################
##
##  LIKELIHOOD
##
################################################################################

blm.likelihood <- function(type = "normal", df = 3)
{
  type <- match.arg(type, c("normal", "t"))
  list(type = type, errorCov = NULL, df = df)
}


blm.prior <- function(fixed.coef = "non-informative", sigma2 = "non-informative")
{
	conjugate <- FALSE
	priorSigma <- sigma2
	priorBeta <- fixed.coef
  
  if(class(priorBeta) != "fbdstn" && priorBeta != "non-informative")
    stop("blm.prior: prior for coefficients is not of the right type")

  else if(class(priorBeta) == "fbdstn") {
    if(!is.element(priorBeta$name, c("normal", "t", "normmix", "tmix")))
      stop("blm.prior: prior specification for coefficients is not supported")
  }

  if(class(priorSigma) != "fbdstn" && priorSigma != "non-informative")
    stop("blm.prior: prior for variance is not of the right type")

  else if(class(priorSigma) == "fbdstn") {
    if(!is.element(priorSigma$name, c("invChisq")))
      stop("blm.prior: prior specification for variance is not supported")
  }

  list(priorBeta = priorBeta, priorSigma = priorSigma, conjugate = conjugate)
}


blm.sampler <- function(nBurnin = 1000, nSamples = 1000, nThin = 1, nChains = 1,
                        init.point = "prior + likelihood", gamma.init = NULL,
                        sigma2.init = NULL)
{
  sigma.init <- sigma2.init
  sampler.type <- "Gibbs"
  number.chains <- nChains
  beta.init <- gamma.init

  if(number.chains <= 0) {
    warning("blm.sampler: zero or negative number of chains requested; using 1 chain")
    number.chains <- 1
  }

  #just in case this is not an integer as it should be  
  if(ceiling(number.chains) - number.chains > 0) {
    new.number.chains <- ceiling(number.chains)
    warning("blm.sampler: non-integer number of chains specified [", number.chains,
            "]; setting number of chains to [", new.number.chains, "]")
    number.chains <- new.number.chains
  }

  beta <- NULL
  sigma <- NULL

  init.point <- match.arg(init.point, c("prior", "prior + likelihood",
                                        "likelihood", "user's choice"))

  if(init.point == "user's choice") {

    if(is.null(beta.init))
      stop("blm.sampler: initial points for beta must be specified")

    if(is.null(sigma.init))
      stop("blm.sampler: initial points for sigma must be specified")

    if(is.vector(beta.init)) {

      if(number.chains > 1) {
        #need the same length (assuming beta dimension = 1)
        if(length(beta.init) != number.chains)
          stop("blm.sampler: there must be ", number.chains, " initial points specified for beta; only ",
                length(beta.init), " points given")
      }
    }

    else if(is.matrix(beta.init)) {
      #need number.chains columns in matrix 
      if(ncol(beta.init)!= number.chains)
        stop("blm.sampler: there must be ", number.chains, " initial points specified for beta")
    }

    if(!is.vector(sigma.init))
      stop("blm.sampler: initial points for sigma must be provided")

    else if(length(sigma.init)!= number.chains)
      stop("blm.sampler: there must be ", number.chains, " initial points specified for sigma")

    beta <- beta.init
    sigma <- sigma.init
  }

  list(control = list(nBurnin = nBurnin, nSamples = nSamples, nThin = nThin),
       number.chains = number.chains,
       init.point = list(beta = beta, sigma = sigma, type = init.point), 
       sampler = sampler.type)
}


