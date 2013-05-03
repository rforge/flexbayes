crosscorr <- function(x, ...)
  UseMethod("crosscorr")


crosscorr.default <- function(x, ...)
  coda::crosscorr(x)


crosscorr.posterior <- function(x, ...)
  coda::crosscorr(x$chains)


################################################################################

effectiveSize <- function(x, ...)
  UseMethod("effectiveSize")


effectiveSize.default <- function(x, ...)
  coda::effectiveSize(x)


effectiveSize.posterior <- function(x, ...)
  coda::effectiveSize(x$chains)


################################################################################

autocorr.plot <- function(x, ...)
  UseMethod("autocorr.plot")


autocorr.plot.default <- function(x, ...)
  coda::autocorr.plot(x, ...)


autocorr.plot.posterior <- function(x, ..., chain = 1)
  coda::autocorr.plot(x$chains[[chain[1]]], ...)


################################################################################

crosscorr.plot <- function(x, ...)
  UseMethod("crosscorr.plot")


crosscorr.plot.default <- function(x, ...)
  coda::crosscorr.plot(x, ...)


crosscorr.plot.posterior <- function(x, ..., chain = 1)
  coda::crosscorr.plot(x$chains[[chain[1]]], ...)


################################################################################

traceplot <- function(x, ...)
  UseMethod("traceplot")


traceplot.default <- function(x, ...)
  coda::traceplot(x, ...)


traceplot.posterior <- function(x, ...)
  coda::traceplot(x$chains, ...)


################################################################################

autocorr <- function(x, ...)
  UseMethod("autocorr")


autocorr.default <- function(x, ...)
  coda::autocorr(x, ...)


autocorr.posterior <- function(x, ...)
  coda::autocorr(x$chains, ...)


################################################################################

densplot <- function(x, ...)
  UseMethod("densplot")


densplot.default <- function(x, ...)
  coda::densplot(x, ...)


densplot.posterior <- function(x, ...)
  coda::densplot(x$chains, ...)


################################################################################

nvar <- function(x)
  UseMethod("nvar")


nvar.default <- function(x)
  coda::nvar(x)


nvar.posterior <- function(x)
  coda::nvar(x$chains)


################################################################################

thin.posterior <- function(x, ...)
{
  x$chains <- coda::thin(x$chains, ...)
  x
}


################################################################################

nchain <- function(x)
  UseMethod("nchain")


nchain.default <- function(x)
  coda::nchain(x)


nchain.posterior <- function(x)
  coda::nchain(x$chains)


################################################################################

niter <- function(x)
  UseMethod("niter")


niter.default <- function(x)
  coda::niter(x)


niter.posterior <- function(x)
  coda::niter(x$chains)


################################################################################

as.mcmc.list.posterior <- function(x, ...)
  as.mcmc.list(x$chains, ...)


################################################################################

gelman.diag <- function(x, ...)
  UseMethod("gelman.diag")


gelman.diag.default <- function(x, ...)
  coda::gelman.diag(x, ...)


gelman.diag.posterior <- function(x, ...)
  coda::gelman.diag(x$chains, ...)


################################################################################

gelman.plot <- function(x, ...)
  UseMethod("gelman.plot")


gelman.plot.default <- function(x, ...)
  coda::gelman.plot(x, ...)


gelman.plot.posterior <- function(x, ...)
  coda::gelman.plot(x$chains, ...)


################################################################################

geweke.diag <- function(x, ...)
  UseMethod("geweke.diag")


geweke.diag.default <- function(x, ...)
  coda::geweke.diag(x, ...)


geweke.diag.posterior <- function(x, ...)
  coda::geweke.diag(x$chains, ...)


################################################################################

geweke.plot <- function(x, ...)
  UseMethod("geweke.plot")


geweke.plot.default <- function(x, ...)
  coda::geweke.plot(x, ...)


geweke.plot.posterior <- function(x, ...)
  coda::geweke.plot(x$chains, ...)


################################################################################

raftery.diag <- function(x, ...)
  UseMethod("raftery.diag")


raftery.diag.default <- function(x, ...)
  coda::raftery.diag(x, ...)


raftery.diag.posterior <- function(x, ...)
  coda::raftery.diag(x$chains, ...)


################################################################################

heidel.diag <- function(x, ...)
  UseMethod("heidel.diag")


heidel.diag.default <- function(x, ...)
  coda::heidel.diag(x, ...)


heidel.diag.posterior <- function(x, ...)
  coda::heidel.diag(x$chains, ...)


################################################################################

cumuplot <- function(x, ...)
  UseMethod("cumuplot")


cumuplot.default <- function(x, ...)
  coda::cumuplot(x, ...)


cumuplot.posterior <- function(x, ...)
  coda::cumuplot(x$chains, ...)


################################################################################


