\name{fbprior}
\alias{fbprior}

\title{FlexBayes Prior Distribution}

\description{
Specify prior distributions for FlexBayes models.
}

\usage{
  fbprior(dstn, ...)
}

\arguments{
  \item{dstn}{a character string naming the distribution.}

  \item{\dots}{parameters for the specified distribution.}
}

\details{

Following is a list of supported prior distributions and their parameters.  The table gives the name of the parameter used in the documentation (taken from Appendix A of Gelman et. al.), the name of the parameter expected by \code{fbprior}, a short description of the type of object expected, and the default value.

\describe{
  \item{normal}{a (multivariate) normal prior distribution for a parameter (vector) characterized by a mean vector \eqn{\mu} and a covariance matrix \eqn{\Sigma}. \tabular{rrlr}{
    Parameter \tab argument \tab description \tab default \cr
    \eqn{\mu} \tab \code{mean} \tab a numeric vector \tab \code{0.0} \cr
    \eqn{\Sigma} \tab \code{Sigma} \tab a positive definite numeric matrix \tab \code{diag} \cr}
  }

  \item{t}{a (multivariate) t prior distribution for a parameter (vector) characterized by a mean vector \eqn{\mu} and a scale matrix \eqn{\Sigma}. \tabular{rrlr}{
    Parameter \tab argument \tab description \tab default \cr
    \eqn{\mu} \tab \code{mean} \tab a numeric vector \tab \code{0.0} \cr
    \eqn{\Sigma} \tab \code{Sigma} \tab a positive definite numeric matrix \tab \code{diag} \cr}
  }

  \item{nonInformative}{a non-informative flat (uniform) prior distribution.}

  \item{beta}{Not yet documented.}

  \item{gamma}{Not yet documented.}

  \item{normal.mixture}{Not yet documented.}

  \item{t.mixture}{Not yet documented.}

  \item{uniformShrinkage}{a non-informative uniform shrinkage prior distribution for a variance or squared scale parameter.  The probability density function is given by \deqn{p(\sigma^{2}) = \sigma_{0}^{2} / (\sigma^{2} + \sigma_{0}^{2})^{2}} \tabular{rrlr}{
    Parameter \tab argument \tab description \tab default \cr
    \eqn{\sigma_{0}} \tab \code{sigma0} \tab a single positive numeric value \tab \code{1.0} \cr}
  }

  \item{duMouchel}{a non-informative \emph{du Mouchel} prior distribution for a variance or squared scale parameter.  The probability density function is given by \deqn{p(\sigma^{2}) = \sigma_{0} / (2 \sigma (\sigma_{0} + \sigma)^{2})} \tabular{rrlr}{
    Parameter \tab argument \tab description \tab default \cr
    \eqn{\sigma_{0}} \tab \code{sigma0} \tab a single positive numeric value \tab \code{1.0} \cr}
  }

  \item{nonInfoPower}{a non-informative power law prior distribution for a variance or squared scale parameter.  The probability density function is given by \deqn{p(\sigma^{2}) \propto \sigma^{2\kappa}} \tabular{rrlr}{
    Parameter \tab argument \tab description \tab default \cr
    \eqn{\kappa} \tab \code{kappa} \tab a single negative numeric value \tab \code{-1.0} \cr}
  }

  \item{massPoint}{a degenerate Dirac delta prior distribution for known variance or squared scale parameters.  \tabular{rrlr}{
    Parameter \tab argument \tab description \tab default \cr
    \eqn{\sigma^{2}} \tab \code{sigma2} \tab a vector of positive numeric values \tab \code{1.0} \cr}
  }

  \item{wishart}{Not yet documented.}

  \item{invWishart}{an inverse Wishart prior distribution for a covariance matrix parameter.  \tabular{rrlr}{
    Parameter \tab argument \tab description \tab default \cr
    \eqn{\nu_{0}} \tab \code{df} \tab a single positive numeric value \tab \code{1.0} \cr
    \eqn{\Sigma} \tab \code{Sigma} \tab a positive definite matrix \tab \code{-1.0} \cr}
  }

  \item{invChisq}{a non-informative power law prior distribution for a variance or squared scale parameter.  The probability density function is given by \deqn{p(\sigma^{2}) \propto \sigma^{-2(\nu_{0}/2 + 1)} \exp{\{-(\nu_{0} \sigma_{0}^{2})/(2 \sigma^{2})\}}} \tabular{rrlr}{
    Parameter \tab argument \tab description \tab default \cr
    \eqn{\nu_{0}} \tab \code{df} \tab a single positive numeric value \tab \code{3} \cr
    \eqn{\sigma_{0}^{2}} \tab \code{sigma0.sq} \tab a single positive numeric value \tab \code{1.0} \cr} Remark: this distribution is the scaled inverse chi-squared distribution described in Gelman, Carlin, Stern, and Rubin (2004).
  }
}
}

\value{
a list containing the following elements:

\item{name}{a character string namimg the distribution.}

\item{parameters}{a list containing the parameters of the distribution.}
}

\section{References}{
Gelman, A., Carlin, J. B., Stern, H.S., and Rubin, D. B. (2004). \emph{Bayesian Data Analysis}, 2nd Edition.  Boca Raton, Florida: Chapman & Hall.

Daniels, M. J. (1999), "A prior for the variance in hierarchical models", \emph{The Canadian Journal of Statistics}, 27: 567-578.

Hobert, J. P., and G. Casella (1996). The effect of improper priors on Gibbs sampling in hierarchical linear mixed models. \emph{J. American Statistical Association}, 91: 1461-1473.
}

\examples{
beta.prior <- fbprior("normal", mean = rep(0, 4), diag(4))
}

\keyword{FlexBayes}
\keyword{prior}


