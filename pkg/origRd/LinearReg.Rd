\name{LinearReg}
\alias{LinearReg}
\title{A linear regression model for validation}

\description{

The linear regression model is  specified as follows. The regression has two
coefficients, namely the intercept beta1 and the slope beta2.  These are given
normal prior distributions, and the variance of the outcome is given an inverse
gamma prior.

The model is as follows, in BUGS format.  The choices of the hyperparameter values are 
arbitrary.  See the \code{CGR} help file for more information.

	{
	
		for( i in 1 : N ) {

			outcome[i] ~ dnorm(mu[i],tau.c)
		
			mu[i] <- beta1 + beta2 * predictor[i]
		
		}
		

		# define the regression coefficients
		
		beta1 ~ dnorm(0, 5)
		
		beta2 ~ dnorm(0, 5)


		# define the outcome precision
		
		tau.c ~ dgamma(3, 1)
		
	}

}
