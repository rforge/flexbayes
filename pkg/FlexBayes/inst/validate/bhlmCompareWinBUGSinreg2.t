{cat("--------------- BUGS vs. bhlm: Linear Regression ----\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhlm
# Description: Compare the results for the bhlm function
#  to those from the BUGS engine, for a linear 
#  regression model fit to the stack data.

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
bugsModel <- function (){
  for( i in 1 : N ) {
    Loss[i] ~ dnorm(mu[i],tau.c)
    mu[i] <- beta[1] + beta[2] * Air.Flow[i] + 
      beta[3]*Water.Temp[i] + beta[4]*Acid.Conc.[i]
  }
  for(j in 1:4){
    beta[j] ~ dnorm(beta.mean, beta.prec)
  }
  tau.c ~ dgamma(a.prec, b.prec)
  sigma <- pow(tau.c, -0.5)

  invisible()		
}

# Create the data frame for fitting the model
stackModelData <- as.list(stack.dat)
stackModelData$N <- length(stack.dat$Loss)

# Set the values for the hyperparameters.  
a.prec <- 1.5
b.prec <- 1.5
beta.mean <- 0.7
beta.prec <- 5.2
stackModelData$a.prec <- a.prec
stackModelData$b.prec <- b.prec
stackModelData$beta.mean <- beta.mean
stackModelData$beta.prec <- beta.prec

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

coefPrior <- bayes.normal(mean.vector=rep(beta.mean,4), 
  covmat = diag(rep(beta.prec^-1, 4)))
  
varPrior <- bayes.invChisq(df=(2*a.prec), 
  sigma0.sq=(b.prec/a.prec))

stackSampler <- bhlm.sampler( nChains=1,
  init.point = "prior", nSamples=1000,
  nThin=5)

stackPrior <- bhlm.prior( fixed.coef = coefPrior,
  error.var = varPrior )
  
bhlmSamples <- bhlm( fixed.formula = Loss ~ ., 
  data = stack.dat, prior = stackPrior,
  sampler = stackSampler)

# Compare the sample distributions.
# Adjust for the multiple test scripts being run.
compareSampleDistributions( 
  getSamples( bhlmSamples ), 
  getSamples( bugsSamples ), print.ks=F, 
  pvalue.cutoff = (0.01 / 100))
}
