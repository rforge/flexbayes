\name{cgrFlexBayes}
\alias{cgrFlexBayes}
\title{Cook-Gelman-Rubin Validation of FlexBayes}

\usage{cgrFlexBayes(modelName = 'LinearReg', verbose = F, 
  engine = 'OpenBUGS', n.CGR.reps = 200, threshold = 0.01)}

\arguments{  
	\item{modelName}{The name of the model to be validated.  The options are 'LinearReg', 'LogisticReg', 
	'LogisticTprior', 'LogisticOverdisperse', 
		'PoissonGLM', 'PoissonTprior', 'PoissonOverdisperse', and 'AEmodel'. Descriptions of some of these 
		models are available by calling \code{help} on the name of the model, for instance \code{help(LinearReg)}.}
	\item{verbose}{Whether to print the p-values for the test.}
	\item{engine}{Which software engine to use for obtaining the posterior samples.  Options are \code{'OpenBUGS'},
	\code{'WinBUGS'}, and \code{'internal'}.  \code{'internal'} refers to the proprietary model-fitting code within 
	the \code{FlexBayes} package.  Some of these options may not be available, depending on the model to be validated.}
	\item{n.CGR.reps}{The number of replications to use in the Cook, Gelman, and Rubin validation test.
	200 is recommended because the sensitivity of the test increases with the number of replications and because 
	more than 200 replications is very computationally intensive.}
	\item{threshold}{If one of the Bonferroni-adjusted p-values is below this threshold, the test fails and a value of 
	FALSE is returned.}
}

\value{
 A \code{T/F} value indicating whether the test passed.  

}

\description{
Validate one of the model-fitting functions in the FlexBayes package, using a validation
technique defined by Cook, Gelman, and Rubin (2006).  For information 
on the validation technique, see the \code{CGR} help file in this package.  
}
