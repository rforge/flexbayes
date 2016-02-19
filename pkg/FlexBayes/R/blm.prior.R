blm.prior <- function(priorBeta = fbprior("nonInformative"),
                      priorSigma = fbprior("nonInformative"))
{
  if(class(priorBeta) != "fbprior")
    stop("blm.prior: prior for coefficients is not of the right type")

  else if(class(priorBeta) == "fbprior") {
    if(!is.element(priorBeta$name, c("nonInformative", "norm", "t", "normmix", "tmix")))
      stop("blm.prior: prior specification for coefficients is not supported")
  }

  if(class(priorSigma) != "fbprior")
    stop("blm.prior: prior for variance is not the right type")

  else if(class(priorSigma) == "fbprior") {
    if(!is.element(priorSigma$name, c("nonInformative", "invChisq")))
    stop("blm.prior: prior specification for variance is not supported")
  }

  list(priorBeta = priorBeta, priorSigma = priorSigma, conjugate = FALSE)
}


