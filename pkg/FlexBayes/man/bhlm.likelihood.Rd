\name{bhlm.likelihood}
\alias{bhlm.likelihood}

\title{Likelihood Specification for a Bayesian Hierarchical Linear Model}

\description{
Set the distribution of the error terms in a Bayesian hierarchical linear model, 
giving rise to the likelihood factor in the posterior distribution.
The error terms can be distributed \bold{normal} or \bold{t}.
}

\usage{
  bhlm.likelihood(type = "normal", df = 3)
}

\arguments{
  \item{type}{a character string. Either \code{"normal"} for normally distributed errors, or \code{"t"} for t distributed errors. error terms are \bold{normally} distributed.}

  \item{df}{a positive integer value specifying the degrees of freedom for t distributed errors.  This argument is ignored when \code{type = "normal"}.}
}

\value{
a list specifying the likelihood structure of a Bayesian hierarchical linear model.
}

\examples{
# default is a normal distribution
my.lklhd <- bhlm.likelihood()

# t-distributed errors
my.lklhd <- bhlm.likelihood(type = "t", df = 3)
}
