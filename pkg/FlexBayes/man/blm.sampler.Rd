\name{blm.sampler}
\alias{blm.sampler}

\title{Specify the Sampler for a Bayesian Linear Model}

\description{Specify the Sampler for a Bayesian Linear Model}

\usage{
blm.sampler(nBurnin = 1000, nSamples = 1000, nThin = 1, nChains = 1,
            start = c("both", "prior", "likelihood", "user"),
            coefs = NULL, sigma = NULL)
}

\arguments{
  \item{nBurnin}{nBurnin}
  \item{nSamples}{nSamples}
  \item{nThin}{nThin}
  \item{nChains}{nChains}
  \item{start}{start}
  \item{coefs}{coefs}
  \item{sigma}{sigma}
}

\keyword{sampler}


