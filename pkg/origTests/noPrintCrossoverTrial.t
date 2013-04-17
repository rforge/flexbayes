{cat("----- No extraneous info printed while fitting crossover trial model --\n");T}
{
# Functions: posteriorSamples
# Description: Make sure that no extraneous info is printed 
# while fitting a crossover trial model
  
  # save printed output to a file
  sink( "temp.txt" )
  
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

  s <- summary(crossoverPost)
  
  sink()
  ( file.lengths("temp.txt")[1] == 0 )
}
