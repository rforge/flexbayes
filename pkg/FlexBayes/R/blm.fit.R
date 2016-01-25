blm.fit <- function(X, Y, errorCov, likDf, sigmaDf, sigmaScale, iType,
                    betaMean, betaCov, betaDf, props, n.mixtures,
                    burnInLength, nSim, sampleFrequency,
                    sampler.type, read.init.point, beta.init, sigma.init,
                    print.stats)
{
  if(n.mixtures == 1)
    props <- 1

  n <- nrow(X)
  p <- ncol(X)
  dim.error.Cov <- length(errorCov)

  #cat("sampler type is ", sampler.type , "  iType is ", iType, "\n")

  keep_tau2_errors <- 0
  #account for beta and sigma
  number.vars <- 2

  if(iType %in% c(2, 5)) #account for tau2 for beta t-prior
    number.vars <- number.vars + 1

  if(iType %in% c(3, 4, 5, 11, 12, 13, 14)) {

    ##t-likelihood case (keep tau2 errors in data augmentation)
    output.samples <- matrix(0.0, nrow = nSim, ncol = p + 1 + n)
    keep_tau2_errors <- 1

    #account for all tau2_errors
    number.vars <- number.vars + n
  }

  else ##exact sampling
  output.samples <- matrix(0.0, nrow = nSim, ncol = p + 1)

  if(print.stats == 1) {

    #the Gibbs sampler case
    if(sampler.type == 0)
    gibbs.stats <- rep(0, number.vars)

    else
    gibbs.stats <- 0

    #the mixture prior case
    if(iType %in% c(7, 8, 9, 10, 11, 12, 13, 14, 15))
    mixture.stats <- rep(0, n.mixtures)

    else
    mixture.stats <- 0
  }

  fit <- .C("fitBayesianLM",
            as.integer(n),
            as.integer(p),
            as.double(Y),
            as.double(t(X)),
            as.integer(dim.error.Cov),
            as.double(errorCov),
            as.double(likDf),
            as.double(sigmaDf),
            as.double(sigmaScale),
            as.integer(iType),
            as.double(betaMean),
            as.double(betaCov),
            as.double(betaDf),
            as.double(props),
            as.integer(n.mixtures),
            as.integer(burnInLength),
            as.integer(nSim),
            as.integer(sampleFrequency),
            as.integer(read.init.point),
            as.double(beta.init),
            as.double(sigma.init),
            as.integer(sampler.type),
            as.integer(print.stats),
            output.samples = as.double(output.samples),
            gibbs.stats = as.double(gibbs.stats),
            mixture.stats = as.double(mixture.stats),
            PACKAGE = "FlexBayes")

  # output is organized as a matrix of nSim x beta dimension
  # the column beyond beta dimension is associated to simulations of sigma
  # (= sqrt(sigma2))
  # if keep tau2 errors then the last n columns contain the simulated
  # sqrt(scale)of each tau2_i

  if(keep_tau2_errors == 0)
    out.samples <- matrix(fit$output.samples, nrow = nSim, ncol = p + 1)
  else
    out.samples <- matrix(fit$output.samples, nrow = nSim, ncol = p + 1 + n)

  list(samples = out.samples,
       gibbs = fit$gibbs.stats,
       mixture = fit$mixture.stats)
}



