profile.blm <- function(fitted, param, ...)
{
  cf.names <- fitted$coef.names
  if(missing(param)) {
    if(length(cf.names) == 1)
      param <- cf.names
    if(length(cf.names) == 2)
      param <- cf.names[2]
  }
  param <- match.arg(param, choices = c(cf.names, "(sigma)"))

  #mle
  mle.m <- fitted$mle$coefficients[param]
  mle.se <- sqrt(fitted$mle$vcov[param, param])
  mle.min <- mle.m - 4.0 * mle.se
  mle.max <- mle.m + 4.0 * mle.se

  #posterior
  posterior.chain <- fitted$chains[[1]][, param]
  posterior.density <- stats::density(posterior.chain)

  #prior
  prior <- fitted$prior$priorBeta$name

  status <- switch(prior,

    "nonInformative" = {
      prior.range <- c(mle.min, mle.max)
      0
    },

    "norm" = {
      prior.m <- fitted$prior$priorBeta$parameters$mean[param]
      prior.sd <- sqrt(fitted$prior$priorBeta$parameters$S[param, param])
      prior.range <- qnorm(c(0.001, 0.999), mean = prior.m, sd = prior.sd)
      0
    },

    "t" = {1},

    "normmix" = {
      prior.m <- sapply(fitted$prior$priorBeta$parameters$mean,
                        function(u, param) u[param],
                        param = param)
      prior.sd <- sapply(fitted$prior$priorBeta$parameters$S,
                         function(u, param) u[param, param],
                         param = param)
      prior.w <- fitted$prior$priorBeta$parameters$w
      prior.nm <- norMix(mu = prior.m, sigma = sqrt(prior.sd), w = prior.w)
      prior.range <- qnorMix(c(0.001, 0.999), prior.nm)
    },

    "tmix" = {1}

  )

  r <- range(c(mle.min, mle.max, prior.range, posterior.density$x))
  x <- seq(from = r[1], to = r[2], length.out = 1000)

  d.prior <- switch(prior,

    "nonInformative" = NA,

    "norm" = dnorm(x, mean = prior.m, sd = prior.sd),

    "t" = {1},

    "normmix" = dnorMix(x, prior.nm),

    "tmix" = {1}
  )

  d.mle <- switch(fitted$mle$distribution,

    "norm" = dnorm(x, mean = mle.m, sd = mle.se),

    "t" = rep(0, 1000)
  )

  upper <- max(c(d.mle, posterior.density$y, d.prior), na.rm = TRUE)

  plot(x, d.mle, type = "l", xlim = r, ylim = c(0.0, upper),
       xlab = param, ylab = "Density", lty = 2)

  lines(posterior.density, lwd = 2)

  if(any(is.na(d.prior))) {
    usr <- par("usr")
    d.prior <- rep(1 / (diff(usr[1:2])), length(x))
  }

  lines(x, d.prior, lty = 3)

  legend("topright", legend = c("Posterior", "Prior", "Likelihood"), lty = c(1, 3, 2),
         lwd = c(2, 1, 1), bty = "n")

  invisible(fitted)
}




