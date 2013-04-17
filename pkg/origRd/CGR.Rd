\name{CGR}
\alias{CGR}
\title{Description of the Cook-Gelman-Rubin validation method}

\source{
S.R. Cook, A. Gelman, D.B. Rubin (2006). Validation of Software for 
Bayesian Models Using Posterior Quantiles.
\emph{Journal of Computational and Graphical Statistics}. 15(3), 675-692.
}

\description{
	One validation method used in the FlexBayes package is that described by 
	Cook, Gelman, and Rubin (2006).  To use this test, a Bayesian model must be specified, including specific 
	values of the hyperparameters.  Then a number of replications of the following procedure are performed.  Parameter
	values are drawn from the prior distribution.  Then data values are drawn
	from the model, given the sampled parameter values.  Posterior samples are then drawn using one of the Bayesian
	inference functions in FlexBayes.  For each parameter the quantile of the true parameter value 
	under the posterior distribution is then estimated using the posterior samples.  
	The authors prove that if the posterior samples are from the correct
	distribution then these quantiles should be uniformly distributed.  
		
	In order
	to check for uniformity, Cook, Gelman, and Rubin (2006) take an inverse normal CDF transformation of the quantiles.  
	The transformed quantiles should then have a normal distribution.  The authors 
	take the sum over the replications of the squares of the transformed quantiles and find a p-value for this
	sum using the right tail of the chi-squared distribution.
	
	The cgrFlexBayes function performs this test, by calling a function written by Cook, Gelman, and Rubin and 
	included by permission in FlexBayes.  
	The cgrFlexBayes function calculates the left-tailed
	p-value as well as the right-tailed p-value.  In addition, 
	it tests the normality of the transformed 
	quantiles directly using a Kolmogorov-Smirinov goodness-of-fit 1-sample
	test.  These three p-values detect different types of departures from uniformity, for the quantiles.  
	
	Using the Cook, Gelman, and Rubin test requires specifying the hyperparameter values for the model.  If hyperparameter values are chosen that imply very diffuse
	priors, then extreme values of the parameters may be sampled from the prior. If this occurs, then the data
	values that are sampled may also be extreme.  For instance, in the Poisson GLM the outcome counts may all be drawn to be zero.  Extreme data values such as these may lead to poor convergence
	of FlexBayes.  To avoid this problem, cgrFlexBayes uses very informative priors
	centered at intermediate values of the parameters
	in all of the models that are tested.  Besides conforming to this requirement, the choice of 
	hyperparameter values in the tested models is arbitrary.  
}
