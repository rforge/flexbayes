\name{LogisticReg}
\alias{LogisticReg}
\title{A logistic regression model for validation}

\description{

The logistic regression model is  specified as follows. The regression has only one
coefficient, namely the intercept beta0.  The risk parameter lambda, instead of being the logistic transformation 
of beta0, is distributed beta with mean equal to the logistic transformation of beta0.  The variance of lambda is 
on the order 1/xi, so that xi is the parameter that adds variability to the logistic regression.  xi is then given
a gamma prior.

The model is as follows, in BUGS format.    The choices of the hyperparameter values are 
arbitrary.  See the \code{CGR} help file for more information.

	{
	

	   for( i in 1 : N ) {
	   
		 	yi[i] ~ dbin (lambda[i], ni[i])
	   
	   	lambda[i] ~ dbeta (alpha[i], beta[i])
			
			alpha[i] <- xi * ( mu[i])
			
			beta[i] <- xi * (1 - mu[i])
		 	
		 	mu[i] <- exp(meanLogitScale[i]) / (1 + exp(meanLogitScale[i]))
			
			meanLogitScale[i] <- beta0	  

	  }


		# the fixed effect coefficient
	   
	   beta0 ~ dnorm( 0, prec.beta)
	   
	   prec.beta ~ dgamma (5000, 50)
	   
	   
		# define the precision parameter 
	   
	   xi ~ dgamma(5000, 50)		
	
	}

}
