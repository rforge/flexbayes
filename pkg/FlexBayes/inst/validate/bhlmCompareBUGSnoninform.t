{cat("------ BUGS vs. bhlm: Linear Regression Improper Priors ----\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhlm
# Description: Compare the results for the bhlm function
#  to those from the BUGS engine, for a linear 
#  regression model fit to the stack data, using 
#  noninformative prior distributions.
#
# The data frame stack.dat combines the stack.loss 
# and stack.x data sets included with S-PLUS.  These 
# data are from the operation of a plant for the 
# oxidation of ammonia to nitric acid, measured on 21 
# consecutive days. 
# The goal is to model  Loss (percent of ammonia lost 
# times 10), as a linear function of Air.Flow (air 
# flow to the plant), Water.Temp (cooling water inlet 
# temperature), and Acid.Conc (acid concentration as 
# a percentage). 

## specify the model in BUGS format	
# Improper priors are not available in BUGS, so 
# approximate using very vague proper priors.
bugsModel <- function (){
  for( i in 1 : N ) {
    Loss[i] ~ dnorm(mu[i],tau.c)
    mu[i] <- beta[1] + beta[2] * Air.Flow[i] + 
      beta[3]*Water.Temp[i] + beta[4]*Acid.Conc.[i]
  }
  for(j in 1:4){
    beta[j] ~ dnorm(0, 1e-6)
  }
  tau.c ~ dgamma(1e-2, 1e-7)
  sigma <- pow(tau.c, -0.5)

  invisible()		
}

# Create the data frame for fitting the model
stackModelData <- stack.dat
stackModelData$N <- length(stack.dat$Loss)

# Specify the parameters for which we wish to save 
# the posterior samples
parameters.to.save <- c("beta", "sigma")

## obtain the posterior samples
bugsSamples <- posteriorSamples (
  data = stackModelData, model = bugsModel, 
  nChains = 1,
  parametersToSave = parameters.to.save,
  nIter = 1000, nBurnin = 1000, nThin = 20,
  engine = "WinBUGS", DIC = F)

# Fit the same model using the bhlm function

stackSampler <- bhlm.sampler( nChains=1,
  init.point = "prior", nSamples=1000,
  nThin=10)

bhlmSamples <- bhlm( fixed.formula = Loss ~ ., 
  data = stack.dat,
  sampler = stackSampler )

# Compare the sample distributions.
# Adjust for the multiple test scripts being run.
compareSampleDistributions( getSamples( bhlmSamples ), 
  getSamples( bugsSamples ), print.ks=F, 
  pvalue.cutoff = (0.01 / 100))
}
