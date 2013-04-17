\name{PoissonGLM}
\alias{PoissonGLM}
\title{A Poisson-outcome GLM for validation}

\source{
M.J. Daniels (1999).  A Prior for the Variance in Hierarchical Models.  \emph{The Canadian Journal
of Statistics} 27(3) 567-578.
}

\description{

The Poisson GLM is specified as follows.   The 
regression has only one coefficient, namely the 
intercept beta0. The risk parameter lambda, instead 
of being the exponential transformation of beta0, is 
distributed gamma with mean equal to the exponential 
transformation of beta0.  The variance of lambda is 
on the order 1/xi, so that xi is the parameter that 
adds variability to the logistic regression.  xi is 
then given a uniform shrinkage prior 
(Daniels, 1999).  The model is as follows, in BUGS 
format.  The choices of the hyperparameter values are 
arbitrary.  See the \code{CGR} help file for more 
information.

	{

	   for( i in 1 : N ) {

	      z[i] ~ dpois(poissonMean[i])

	      poissonMean[i] <- lambda[i] * e[i]

	      lambda[i] ~ dgamma(xi,invScale)

	   }

		# Define the precision parameter for the risk
		
	   xi <- (100 * shrink) / (1 - shrink)

	   shrink ~ dunif(0,1)

	  # Define the mean parameter for the risk
	  
	   mu <- exp(beta0)

	   invScale <- xi / mu

	   beta0 ~ dnorm( 0,prec)

	   prec ~ dgamma( 2.5, 0.25)
	   
	}

}

