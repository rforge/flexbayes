\name{AEmodel}
\alias{AEmodel}
\title{The adverse event model}

\source{
Berry, SM and Berry, DA. (2004) Accounting for Multiplicities in Assessing
Drug Safety. Biometrics. 60. p418-426.
}

\description{

The Berry and Berry (2004) adverse event model is specified as follows, in BUGS format.  The choices of the hyperparameter values are 
arbitrary.  See the \code{CGR} help file for more information.

{

  for (i in 1:B) {

    for(j in 1:Nobs[i]) {
      
      X[i, j] ~ dbin(c[i, j], Nc)

      Y[i, j] ~ dbin(t[i, j], Nt)
      
      logit(c[i,j] )<- gamma[i,j]
      
      logit(t[i,j] )<- gamma[i,j] + theta[i,j]
      
      gamma[i,j]~dnorm(mu.gamma[i], sigma2inv.gamma) 
      
      theta[i,j] <- (1-noDrugEffect[i,j]) * theta.latent[i,j] 
      
      noDrugEffect[i,j] ~ dbern(pi[i])
      
      theta.latent[i,j] ~ dnorm(mu.theta[i], sigma2inv.theta[i])

    } 

    mu.gamma[i] ~ dnorm( mu.gamma.0, tau2inv.gamma)
    
    mu.theta[i] ~ dnorm( mu.theta.0, tau2inv.theta)
    
    pi[i] ~ dbeta (alpha.pi, beta.pi)
    
    sigma2inv.theta[i] ~ dgamma(alpha.sigma.theta, beta.sigma.theta)
    
    sigma.theta[i] <- 1/sqrt(sigma2inv.theta[i])

  }
		
  ## The prior for sigma2inv.gamma
			
  sigma2inv.gamma ~ dgamma(alpha.sigma.gamma, beta.sigma.gamma)   
  
  sigma.gamma <- 1/sqrt(sigma2inv.gamma)
  
  ## The priors for tau2inv.theta and tau2inv.gamma
  
  tau2inv.theta ~ dgamma(alpha.tau.theta, beta.tau.theta)
  
  tau2inv.gamma ~ dgamma(alpha.tau.gamma, beta.tau.gamma)
  
  tau.gamma <- 1/sqrt(tau2inv.gamma)
  
  tau.theta <- 1/sqrt(tau2inv.theta)
  
  mu.gamma.0 ~ dnorm(mu.gamma.00, tau2inv.gamma.00)
  
  mu.theta.0 ~ dnorm(mu.theta.00, tau2inv.theta.00)
		
  ## The priors for alpha.pi and beta.pi
		
  alpha.pi <- alphaExp + 1

  alphaExp ~ dexp(lambda.alpha)

  beta.pi <- betaExp + 1

  betaExp ~ dexp(lambda.beta)

}

}
