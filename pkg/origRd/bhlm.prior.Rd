\name{bhlm.prior}
\alias{bhlm.prior}
\title{
  Prior Specification for Bayesian Hierarchical Linear Models
}
\description{
Specify the prior distributions
for the coefficient and variance parameters of a Bayesian hierarchical linear model, to be fit
through a call to \code{\link{bhlm}}. The model can be written as follows:

  Y_ij = X_ij beta_j + M_ij gamma + e_ij

  beta_j = Z_j alpha + u_j   \bold{:the random effects}

  j = 1, ..., J: the group index
    
  i = 1, ..., n_j: the observation index within group j  
    
  e_ij ~ p( e_ij; sigma^2_j )   \bold{:the error distribution.  May be normal or t.}

  u_jk ~ N( 0, tau^2_k )  \bold{OR:}  u_j ~ N( 0, V )  

  gamma ~ p( gamma )   \bold{:the fixed effects.}

  alpha ~ p( alpha )   \bold{:the second-level effects.}

  where the tau^2_k are given univariate priors or V is given an inverse Wishart prior.
  
  The distribution p( gamma ) can be equal to N( gamma_0, V1_0 ), t( df_gamma, gamma_0, V1_0 ), or 
  it can be flat (improper).  The same holds for p( alpha ).
}
\usage{
bhlm.prior(error.var=bayes.nonInfoPower(-1.), 
  fixed.coef="non-informative", 
  level2.coef="non-informative", 
  random.var=bayes.nonInfoPower(-1.), 
  common.error.var=2)
}
\arguments{
\item{error.var }{
an object of class \code{bayes.distribution} specifying
the prior distribution for the error variance \bold{sigma^2}.
The only available informative prior distribution for \bold{sigma^2} is
 inverse chi-squared, specified via a call to
 \code{\link{bayes.invChisq}}. One can also specify a fixed known value for \bold{sigma^2}
 via a call to \code{\link{bayes.massPoint}}.

Available non-informative priors for \bold{sigma^2} are: power law (a negative power of \bold{sigma^2}),
uniform shrinkage, and DuMouchel, specified via a call to 
\code{\link{bayes.nonInfoPower}}, 
\code{\link{bayes.uniformShrinkage}}, 
or \code{\link{bayes.duMouchel}}, respectively.

If \code{common.error.var} = 2 (the default), then all groups are assumed to have the same
error variance \bold{sigma^2}. If \code{common.error.var} = 1,
 then the groups have different error variances, having a common prior distribution. 
If \code{common.error.var} = 0, then the groups have different error variances with different priors. In the
latter case \code{error.var} must be a list of number-of-groups objects of
 class \code{bayes.distribution}, where
each element of the list specifies the prior for \bold{sigma^2_j}. All these priors must
be of the same type; only the (hyper)parameters of the prior may differ.
}
\item{fixed.coef }{
either the string \code{"non-informative"}, 
or an object of class \code{bayes.distribution} specifying
the prior distribution for the fixed effect coefficients (\bold{gamma}).
The default distribution is a non-proper flat prior specified via the string 
\code{"non-informative"}. Available informative priors for
 the fixed effect coefficients are
 \bold{normal} or \bold{t}, specified via a call to 
\code{\link{bayes.normal}} or
\code{\link{bayes.t}}, respectively. 
}
\item{level2.coef }{
either the string \code{"non-informative"}, 
or an object of class \code{bayes.distribution} specifying
the prior distribution for the second level effect coefficients (\bold{alpha}).
The default distribution is an improper flat prior specified via the string 
\code{"non-informative"}. Available informative priors for
 the second-level effect coefficients are
\bold{normal} or \bold{t}, specified via a call to 
\code{\link{bayes.normal}} or
\code{\link{bayes.t}}, respectively. 
}
\item{random.var }{
an object of class \code{bayes.distribution} specifying
the prior distribution for the variance / covariance structure
 of the random effect coefficients.
Either the random effects \bold{beta_jk} are assumed to be independently distributed 
according to normal distributions with variance \bold{tau^2_k}, or the random effect vectors \bold{beta_j} are assumed
to be independently distributed according to multivariate normal distributions with 
covariance matrix \bold{V}.  
Accordingly, one can specify \code{random.var} as an object of class \code{bayes.distribution}, specifying a common 
univariate prior for \bold{tau^2_k}, as a list of objects of class \code{bayes.distribution}, specifying the 
univariate priors for each \bold{tau^2_k}, or as  an object of class \code{bayes.distribution} that specifies a
prior for \bold{V}.
The only univariate informative prior distribution for \bold{tau^2_k} is inverse chi-squared, specified via a call to
\code{\link{bayes.invChisq}}. 
Available non-informative univariate priors for \bold{tau^2_k} are power law (the default),
uniform shrinkage, and DuMouchel, specified via a call to
\code{\link{bayes.nonInfoPower}}, 
\code{\link{bayes.uniformShrinkage}}, 
or \code{\link{bayes.duMouchel}}, respectively.
The only available prior for \bold{V} is inverse Wishart, specified via a call to 
\code{\link{bayes.invWishart}}. 
}
\item{common.error.var }{
an integer specifying the structure of the error variances \bold{sigma^2_j} for the groups.
Possible values of \bold{ common.error.var } are:

 \bold{2:} all groups are assumed to have the same
error variance \bold{sigma^2}. 

 \bold{1:} the groups have different error variances, having a common prior distribution. 

 \bold{0:} the groups have different error variances with different priors. In this
 case \code{error.var} 
must be a list of number-of-groups objects of
 class \code{bayes.distribution}, where each element of the list specifies the prior distribution for \bold{sigma^2_j}. All priors must
be of the same type; only the (hyper)parameters of the prior may differ.
}
}
\value{
a list specifying the prior distributions for a Bayesian hierarchical linear model.
}
\section{References}{

Gelman, A., J. B. Carlin, H. S. Stern, and D. B. Rubin (2004). \emph{ Bayesian Data Analysis}, 2nd Edition. Boca Raton, Florida: Chapman & Hall.

Morris, C. N., and S. L. Normand, "Hierarchical models for combining information and for meta-analyses", In \emph{Bayesian Statistics 4, J. M. Bernardo, J. O. Berger, A. P. Dawid, and F. M. Smith Eds.},  pages 321-344, Oxford University Press, 1992.

Hobert, J. P., and G. Casella, "The effect of improper priors on Gibbs sampling in hierarchical linear mixed models", \emph{J. American Statistical Association}, vol. 91, 436: 1461- 1473, 1996.

Daniels, M. J., "A prior for the variance in hierarchical models", \emph{The Canadian Journal of Statistics}, vol. 27, 3: 567-578, 1999.

Wakefield, J. C., F. M. Smith, A. Racine-Poon, and A. E. Gelfand, "Bayesian analysis of linear and non-linear population models by using the Gibbs sampler", \emph{Journal of the Royal Statistical Society Series C (Appl. Statist.)}, vol. 43, 1: 201-221, 1994.
}
\seealso{
\code{\link{bhlm}},
\code{\link{bhlm.likelihood}},
\code{\link{bhlm.sampler}}.
}
\examples{

alpha.prior <- bayes.nonInformative()
error.prior <- bayes.invChisq(3,1)
betaCov.prior <- bayes.invChisq(3,1)

orthoPrior <- bhlm.prior(error.var=error.prior, 
  level2.coef = alpha.prior, 
  random.var=betaCov.prior, common.error.var=2)

orthoLike <- bhlm.likelihood(type="normal")

orthoSampler <- bhlm.sampler( nBurnin=1000, 
  nSamples=1000, nThin=50,
  nChains=1, init.point = "prior")

bhlm.samples <- bhlm(random.formula = distance~I(age-11),
  level2.formula = ~Sex, group.formula = ~Subject, 
  data = Orthodont, prior = orthoPrior, 
  likelihood = orthoLike, sampler = orthoSampler)

}
\keyword{Bayes}
% docclass is function
% Converted by Sd2Rd version 37351.
