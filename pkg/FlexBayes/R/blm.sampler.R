blm.sampler <- function(nBurnin = 1000, nSamples = 1000, nThin = 1, nChains = 1,
                        start = c("both", "prior", "likelihood", "user"),
                        coefs = NULL, sigma = NULL)
{
  start <- match.arg(start)
  nChains <- round(nChains)

  if(nChains <= 0) {
    warning("zero or negative number of chains requested; using 1 chain")
    nChains <- 1
  }

  if(start == "user") {
    if(is.null(coefs))
      stop("no initial point for coefficients specified")

    if(is.null(sigma))
      stop("initial points for sigma must be specified")

    if(is.vector(coefs)) {
      if(nChains > 1) {
        #need the same length (assuming beta dimension = 1)
        if(length(coefs) != nChains)
          stop("there must be ", nChains, " initial points specified ",
               "for the coefficients; only ", length(coefs),
               " points given")
      }
    }

    else if(is.matrix(coefs)) {
      #need nChains columns in matrix
      if(ncol(coefs)!= nChains)
        stop("there must be ", nChains, " initial points specified for ",
             "the coefficients")
    }

    if(!is.vector(sigma.init))
      stop("initial points for sigma must be provided")

    else if(length(sigma.init)!= nChains)
      stop("there must be ", nChains, " initial points specified for sigma")
  }

  list(control = list(nBurnin = nBurnin, nSamples = nSamples, nThin = nThin),
       nChains = nChains,
       init.point = list(coefs = coefs, sigma = sigma, type = start),
       sampler = "Gibbs")
}
