\name{bhlm.likelihood}
\alias{bhlm.likelihood}
\title{
  Likelihood Specification for a Bayesian Hierarchical Linear Model
}
\description{
Set the distribution of the error terms in a Bayesian hierarchical linear model, 
giving rise to the likelihood factor in the posterior distribution.
The error terms can be \bold{normally} or \bold{t} distributed. 
}
\usage{
bhlm.likelihood(type="normal", df=3)
}
\arguments{
\item{type}{
a string specifying the likelihood type.  Options are:

 \bold{"normal":} error terms are \bold{normally} distributed. 

 \bold{"t":} error terms are \bold{t}-distributed. 
}
\item{df}{
an integer specifying the \bold{degrees of freedom} associated to the
 \bold{t}-distributed errors.  If \code{type = "normal"}, this argument
 is ignored.
}
}
\value{
a list specifying the likelihood structure of a Bayesian hierarchical linear model, to be passed to a call to \code{bhlm}.
}
\seealso{
\code{\link{bhlm.sampler}}.
}
\examples{
# default is a normal distribution
my.lklhd <- bhlm.likelihood()

# t-distributed errors
my.lklhd <- bhlm.likelihood( type = "t", df = 3 )

# Fit the Orthodont data
alpha.prior <- bayes.normal(mean=zero, cov=identity)
error.prior <- bayes.invChisq(3,1)
betaCov.prior <- bayes.invWishart( df=3,scale=diag(2) )

orthoPrior <- bhlm.prior(error.var=error.prior, 
  level2.coef = alpha.prior, 
  random.var=betaCov.prior, common.error.var=2)

orthoLike <- bhlm.likelihood(type = "normal")

orthoSampler <- bhlm.sampler( nBurnin=10000, 
  nSamples=1000, nThin=50,
  nChains=1, init.point = "prior")

bhlm.samples <- bhlm(random.formula = distance~I(age-11),
  level2.formula = ~Sex, group.formula = ~Subject, 
  data = Orthodont, prior = orthoPrior, 
  likelihood = orthoLike, sampler = orthoSampler)

effectiveSize(bhlm.samples)
summary(bhlm.samples, maxVars=60)
}
