\name{runAEModel}
\alias{runAEModel}
\title{Run a Model for Adverse Events}

\usage{
runAEModel(aeCounts, nIter = 5000, nBurnin = 10000,
	nThin = 10, diffusePrior = F, piSymm = T,
	engine = "OpenBUGS")
}
\arguments{
\item{aeCounts}{The adverse events count data set.  This should be a list
with the following elements:

B: The number of body systems

Nt: The number of patients in the treatment group

Nc: The number of patients in the control group

Nobs: A vector of the number of adverse event types in each body system.

Y: A matrix with B rows and max(Nobs) columns.  Each entry (i,j) is the number
of adverse event occurrences in the treatment group for the jth AE type in the ith body system.  If
Nobs[i] < j, then entry (i,j) should be equal to NA.  Multiple occurrences of the adverse event in
a particular patient is counted only once, so that entry (i,j) should always be less than Nt.

X: A matrix with B rows and max(Nobs) columns.  Each entry (i,j) is the number
of adverse event occurrences in the control group for the jth AE type in the ith body system.  If
Nobs[i] < j, then entry (i,j) should be equal to NA.  Multiple occurrences of the adverse event in
a particular patient is counted only once, so that entry (i,j) should always be less than Nc.}
\item{nIter}{The number of iterations of the Markov chain sampler, after burnin and thinning.}
\item{nBurnin}{length of burn in, i.e. number of iterations to discard at the beginning.}
\item{nThin}{thinning rate.  Must be a positive integer.}
\item{diffusePrior}{Whether to use a more diffuse version of the prior from Berry and Berry (2004) (if TRUE)
or to use the original hyperparameter values from Berry and Berry (2004) (if FALSE).  Results of model-fitting should
be checked to verify that they do not change dramatically using these different prior choices.}
\item{piSymm}{The choices of the hyperparameter values in Berry and Berry (2004) yield a symmetric prior 
distribution on the parameter pi_b as defined in that paper.  This implies that the prior probability of
a treatment effect for body system b is 0.5.  This may not be desirable in the sense that for many 
pharmaceuticals, a priori we do not expect that they affect half of the body systems.  piSymm = F sets a much 
lower prior probability of a treatment effect for each body system.  Model results should be checked to verify that 
they are robust to this choice.}
\item{engine}{which software engine to use to generate the posterior samples.}
}
\description{
	Run the model for clinical adverse events given in Berry and Berry (2004).
}
\source{
Berry, SM and Berry, DA. (2004) Accounting for Multiplicities in Assessing
Drug Safety. Biometrics. 60. p418-426.
}
\details{
For an example, see the \code{drugAdverseEvents} help file. 
}

\value{
	An object of class "posterior" that contains the posterior samples from the Berry and Berry (2004) model.  
}

