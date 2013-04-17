\name{bhbm.prior}
\alias{bhbm.prior}
\title{
  Prior Specification for Bayesian Hierarchical Binomial Models
}
\description{
Set the prior distribution specification 
for the coefficient and variance parameters of a Bayesian hierarchical binomial model,
to be fit through a call to \code{\link{bhbm}}. The models can be written as follows:

\code{    for j = 1, ..., J: the group index}

\code{    and i = 1, ..., n_j: the index for observations within group j,}

 \bold{No overdispersion case.}

\code{    counts_ij ~ Binomial( n_ij, theta_ij )  }

\code{    logit( theta_ij ) = X_ij beta_j + M_ij gamma  }\bold{:the link function}

\code{    beta_j = Z_j alpha + u_j   } \bold{:the random effects}

\code{    gamma ~ p( gamma )} \bold{:the fixed effects}

\code{    alpha ~ p( alpha )} \bold{:the second-level effects}

\code{  u_jk ~ N( 0, tau^2_k )  \bold{OR:}  u_j ~ N( 0, V )  }

\bold{ where}

\code{    counts_ij: data response (count) for the i-th case in the j-th group}

\code{    n_ij: number of trials for the i-th case in the j-th group}

\code{    X_ij: data predictors for the random effects of the i-th case in the j-th group}

\code{    M_ij: data predictors for the fixed effects of the i-th case in the j-th group}

\code{    Z_j: data predictors for the "population" or second-level effects in the j-th group}

\code{  The tau^2_k are given univariate priors or V is given an inverse Wishart prior.}

\code{    The distribution p( gamma ) can be equal to N( gamma_0, V1_0 ) or t( df, gamma_0, V1_0 ), or 
  it can be flat (improper).  The same holds for p( alpha ). }


 \bold{Beta conjugate overdispersion.}

\code{    counts_ij ~ Binomial( n_ij, theta_ij )  }

\code{    theta_ij ~ Beta( xi_j mu_ij, xi_j ( 1 -  mu_ij ) ) }

\code{    logit( mu_ij ) = X_ij beta_j + M_ij gamma  }\bold{:the link function}

\code{    xi_j ~ Uniform Shrinkage( xi_0 )}

    and the rest of the hierarchy is as given above.


\bold{Logit-normal overdispersion.}

\code{     counts_ij ~ Binomial( n_ij, theta_ij )  }

\code{     logit( theta_ij )~ Normal(  mu_ij, sigma2_j ) }

\code{     mu_ij = X_ij beta_j + M_ij gamma   } \bold{:the link function }

\code{     sigma2_j ~ p( sigma2_j ): can be common for all groups; or known, etc. }

   and the rest of the hierarchy is as given above.


}
\usage{
bhbm.prior( sigma2=NULL, xi = NULL, 
  fixed.coef="non-informative", 
  level2.coef="non-informative", 
  random.var=bayes.nonInfoPower( -1.0 ), 
  common.glm=3)
}
\arguments{
\item{sigma2}{
an object of class \code{bayes.distribution} specifying
the prior distribution for the parameter \bold{sigma^2} (in the binomial logit-normal model).  
The only acceptable informative prior distribution for \bold{sigma^2} is the 
inverse Chi-squared distribution, specified through 
\code{\link{bayes.invChisq}}.
One can also specify a fixed known value for \bold{sigma^2}
 through \code{\link{bayes.massPoint}}.

Acceptable non-informative priors for \bold{sigma^2} are non-informative power,
uniform shrinkage, and du Mouchel, specified through 
\code{\link{bayes.nonInfoPower}}, 
\code{\link{bayes.uniformShrinkage}}, 
and \code{\link{bayes.duMouchel}}, respectively.

If \code{common.glm = 2}, then all groups are assumed to have the same
variance \bold{sigma^2}. If \code{common.glm = 1},
 then the groups are not
assumed to have the same variance, but they share the prior for it. 
If \code{common.glm = 0}, then the groups are not
assumed to have the same variance, nor share the same prior for it. Hence, in this
latter case,  \code{sigma2} must be a list of [number of
 groups] objects of
 class \code{bayes.distribution}.
Each element of the list specifies the prior distribution for \bold{sigma^2_j}. All priors must
be of the same type, only the hyperparameters in the prior may differ.
}

\item{xi}{
an object of class \code{bayes.distribution} specifying
the prior distribution for the parameter \bold{xi} (in 
the binomial beta-conjugate model); the only valid prior distribution is uniform shrinkage,
specified through a call to \code{\link{bayes.uniformShrinkage}}.

If \code{common.glm = 2}, then all groups are assumed to have the same
precision \bold{xi}. If \code{common.glm = 1},
 then the groups are not
assumed to have the same precision, but they share the prior for it. 
If \code{common.glm = 0}, then the groups are not
assumed to have the same precision, nor share the same prior for it. Hence, in this
latter case,  \code{xi} must be a list of [number of
 groups] objects of
 class \code{bayes.distribution}.
Each element of the list specifies the prior distribution for \bold{xi_j}.
}

\item{fixed.coef }{
either the string \code{"non-informative"}, 
or an object of class \code{bayes.distribution} specifying
the prior distribution for the fixed effects coefficients \bold{gamma}.
The default distribution is an improper flat prior specified through the string 
\code{"non-informative"}. Acceptable informative priors for
 the fixed effects coefficients are
\bold{normal} or \bold{t}, specified through 
\code{\link{bayes.normal}} and
\code{\link{bayes.t}}, respectively. 
}
\item{level2.coef }{
either the string \code{"non-informative"}, 
or an object of class \code{bayes.distribution} specifying
the prior distribution for the second level effects coefficients (\bold{alpha}) of a Bayesian 
hierarchical model.
The default distribution is an improper flat prior specified through the string 
\code{"non-informative"}. Available informative priors for
 the fixed coefficients are
\bold{normal} or \bold{t}, specified through 
\code{\link{bayes.normal}} and
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
\item{common.glm }{
an integer specifying whether overdispersion will be used, and if so 
whether the overdispersion parameter \bold{sigma^2} or \bold{xi} is group-specific and the type of 
prior to be used for this parameter.

Acceptable values of \code{ common.glm } are:

 \bold{3:} no overdispersion. 

 \bold{2:} all groups are assumed to have the same
overdispersion parameter \bold{sigma^2} or \bold{xi}. 

 \bold{1:} the groups are not
assumed to have the same overdispersion parameter, but they share the prior for it. 

 \bold{0:} the groups are not
assumed to have the same overdispersion, nor share the same prior for it. Hence, in this
latter case, the arguments \code{xi} or \code{sigma2} 
must be a list of [number of groups] objects of
 class \code{bayes.distribution}.
Each element of the list specifies the prior distribution for \bold{sigma^2_j} or \bold{xi_j}. All priors must
be of the same type, only the hyperparameters in the prior may differ.
}
}
\value{
a list specifying the prior information of a Bayesian hierarchical binomial model.
}
\section{References}{

Gelman, A., J. B. Carlin, H. S. Stern, and D. B. Rubin (2004). \emph{ Bayesian Data Analysis}, 2nd Edition. Boca Raton, Florida: Chapman & Hall.

Christiansen, C. L., and C. N. Morris, "Hierarchical Poisson regression modeling", \emph{J. American Statistical Association}, vol. 92, 438: 618-632, 1997.

Kass, R. E., and D. Steffey, "Approximate Bayesian inference in conditionally independent hierarchical models (parametric empirical Bayes models)", \emph{J. American Statistical Association}, vol. 84, 407: 717-726, 1989.
}
\seealso{
\code{\link{bhpm}},
\code{\link{bhbm}},
\code{\link{bhbm.sampler}}.
}
\examples{
my.prior <- bhbm.prior( xi = bayes.invChisq( df = 3, sigma0.sq = 1 ),
  common.glm = 2 )
#or 
my.prior <- bhbm.prior( sigma2 = bayes.uniformShrinkage( 2 ),
  common.glm = 1 )

}
\keyword{Bayes}
