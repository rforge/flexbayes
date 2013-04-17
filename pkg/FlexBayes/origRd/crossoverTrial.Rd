\name{crossoverTrial}
\alias{crossoverTrial}
\docType{data}
\title{Crossover drug trial data}
\usage{data(crossoverTrial)}

\description{
  Data from a two-period crossover trial, as given in Gelfand et al. (1990).
  The trial was designed to evaluate the relative efficacy of two drugs.  
  Here we fit a model similar to that from Gelfand et al. (1990), assess convergence, and obtain
  posterior inferences.  Our modification of their model uses prior distributions for the variance 
  components that are uniform on a large interval for the standard deviation parameter, as suggested
  in Gelman, Carlin, Stern, and Rubin (2004).
}

\examples{

# Here we fit a hierarchical linear model using the 
# posteriorSamples function, although the same model
# can also be fit using the bhlm function (with 
# slightly different prior distributions).

# Specify the model in BUGS format
bugsModel <- function(){
  for (i in 1:N){
    y[i] ~ dnorm (mean[i], prec1)
    mean[i] <- 
      mu + patientEff[i] + formEff[i] + periodEff[i]
    patientEff[i] <- delta[patient[i]]
    formEff[i] <- phi/2* pow(-1, form[i])
    periodEff[i] <- pi/2* pow(-1, period[i])
  }
  for (j in 1:nPatients){
    delta[j] ~ dnorm(0, prec2)
  }
  mu ~ dnorm(0, 0.01)
  phi ~ dnorm(0, 0.01)
  pi ~ dnorm(0, 0.01)
  prec1 <- pow(sigma1, -2)
  prec2 <- pow(sigma2, -2)
  sigma1 ~ dunif(0, 100)
  sigma2 ~ dunif(0, 100)
  theta <- exp(phi)
  
  invisible()
}

# Get the data
data(crossoverTrial)

# Specify the initial values
drugInits <- list(mu=0, phi=0, pi=0, sigma1= 1, 
  sigma2 = 1)
drugInits <- list(drugInits)

# Specify the parameters to save	
paramsToSave <- c ("mu", "phi", "pi", "prec1", 
  "prec2", "theta")

# Obtain the posterior samples.  THIS TAKES A FEW
# MINUTES
crossoverPost <- 
  posteriorSamples(model = bugsModel, 
    data = crossoverTrial, 
    inits = drugInits, parametersToSave = paramsToSave, 
    nIter = 10000, nThin = 20, engine = "WinBUGS")

traceplot(crossoverPost)		
# The traceplots indicate convergence.

autocorr.plot(crossoverPost)
# Here, enough thinning has been done
# so that the resulting autocorrelations are negligible.

geweke.plot(crossoverPost)
# The Geweke-Brooks plots are acceptable.

# since the convergence diagnostics look good, obtain 
# the posterior estimates and intervals.
summary(crossoverPost)
densplot(crossoverPost)

# The posterior 95-percent interval for the drug effect 
# contains zero.  Therefore there is no strong 
# evidence of a drug effect.  However, the period 
# effect interval is entirely below zero, so that 
# there is strong evidence for a period effect.
}

\source{
  A. E. Gelfand, S. E. Hills, A. Racine-Poon, and A. F. M. Smith (1990):
  Illustration of Bayesian Inference in Normal Data Models Using Gibbs Sampling,
  \emph{Journal of the American Statistical Association} 85 (412), 972-985.
  
  Gelman, A., Carlin, J.B., Stern, H.S., Rubin, D.B. (2004): 
  \emph{Bayesian Data Analysis}, 2nd edition, CRC Press.

}
