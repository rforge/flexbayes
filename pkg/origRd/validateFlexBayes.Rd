\name{validateFlexBayes}
\alias{validateFlexBayes}
\title{Validate flexBayes}

\usage{validateFlexBayes(modelName = 'LinearReg', verbose = F, 
  engine = 'OpenBUGS', n.CGR.reps = 200, threshold = 0.01)}

\arguments{  
	\item{modelName}{The name of the model to be validated.  The options are 'LinearReg', 'LogisticReg', 
	'LogisticTprior', 'LogisticOverdisperse', 
		'PoissonGLM', 'PoissonTprior', 'PoissonOverdisperse', and 'AEmodel'. Descriptions of some of these 
		models are available by calling \code{help} on the name of the model, for instance \code{help(LinearReg)}.}
	\item{verbose}{Whether to print the p-values for the test or simply return a \code{T} / \code{F} result.}
	\item{engine}{Which software engine to use for obtaining the posterior samples.  Options are \code{'OpenBUGS'},
	\code{'WinBUGS'}, and \code{'internal'}.  \code{'internal'} refers to the proprietary model-fitting code within 
	the \code{FlexBayes} package.  Some of these options may not be available, depending on the model to be validated.}
	\item{n.CGR.reps}{The number of replications to use in the Cook, Gelman, and Rubin validation test.
	200 is recommended because the sensitivity of the test increases with the number of replications and because 
	more than 200 replications is very computationally intensive.}
	\item{threshold}{If one of the Bonferroni-adjusted p-values is below this threshold, the test fails.  Only
	used if \code{verbose = F}.}
}

\value{
 A \code{T/F} value indicating whether the test passed.  If a model is requested that has not been defined, a \code{NULL} value is returned.  

}

\description{
Validate one of the model-fitting functions in the FlexBayes package.  The validation technique is 
a modification of that defined by Cook, Gelman, and Rubin (2006).  For information 
on the validation technique, see the \code{CGR} help file in this package.  

If \code{verbose = T}, the p-values for several tests are printed.
}
