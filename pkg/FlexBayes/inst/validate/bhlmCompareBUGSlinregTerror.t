{cat("----- BUGS vs. bhlm: Linear Regression T Errors ----\n");T}
{
# Functions: posteriorSamples(engine='WinBUGS'), bhlm
# Description: Compare the results for the bhlm function
#  to those from the BUGS engine, for a linear 
#  regression model with t errors fit to the stack data.

# The data frame stack.dat combines the stack.loss 
# and stack.x data sets included with S-PLUS.  These 
# data are from the operation of a plant for the 
# oxidation of ammonia to nitric acid, measured on 21 
# consecutive days. 
# The goal is to model Loss (percent of ammonia lost 
# times 10), as a linear function of Air.Flow (air 
# flow to the plant), Water.Temp (cooling water inlet 
# temperature), and Acid.Conc (acid concentration as 
# a percentage). 

## specify the model in BUGS format	
# specify the t distributed errors using a scale 
# mixture of normals.
bugsModel <- function (){
  for( i in 1 : N ) {
    Loss[i] ~ dnorm(mu[i], prec[i])
    prec[i] <- extradisperse[i] * sigma2inv / df.Loss
    extradisperse[i] ~ dchisqr(df.Loss)
    tau[i] <- sqrt(df.Loss / extradisperse[i])
    mu[i] <- beta[1] + beta[2] * Air.Flow[i] + 
      beta[3]*Water.Temp[i] + beta[4]*Acid.Conc.[i]
  }
  for(j in 1:4){
    beta[j] ~ dnorm(beta.mean, beta.prec)
  }
  sigma2inv ~ dgamma(a.sigma2inv, b.sigma2inv)
  sigma <- pow(sigma2inv, -0.5)

  invisible()		
}

# Create the data frame for fitting the model
stackModelData <- stack.dat
stackModelData$N <- length(stack.dat$Loss)

# Set the values for the hyperparameters.  
a.sigma2inv <- 0.6
b.sigma2inv <- 4.2
beta.mean <- 0.7
beta.prec <- 0.15
df.Loss <- 4.1
stackModelData$a.sigma2inv <- a.sigma2inv
stackModelData$b.sigma2inv <- b.sigma2inv
stackModelData$beta.mean <- beta.mean
stackModelData$beta.prec <- beta.prec
stackModelData$df.Loss <- df.Loss

# Specify the parameters for which we wish to save 
# the posterior samples
parameters.to.save <- c("beta", "sigma", "tau")

## obtain the posterior samples
bugsSamples <- posteriorSamples (
  data = stackModelData, model = bugsModel, 
  nChains = 1,
  parametersToSave = parameters.to.save,
  nIter = 1000, nBurnin = 10000, nThin = 500,
  engine = "WinBUGS", DIC = F)
  
# Fit the same model using the bhlm function

coefPrior <- bayes.normal(mean.vector=rep(beta.mean,4), 
  covmat = diag(rep(beta.prec^-1, 4)))
  
varPrior <- bayes.invChisq(df=(2*a.sigma2inv), 
  sigma0.sq=(b.sigma2inv/a.sigma2inv))

stackSampler <- bhlm.sampler( nChains=1, nThin=5,
  init.point = "prior", nSamples=1000 )

stackLikelihood <- bhlm.likelihood(type = "t", df=df.Loss)

stackPrior <- bhlm.prior(fixed.coef = coefPrior, 
  error.var = varPrior)

bhlmSamples <- bhlm( fixed.formula = Loss ~ ., 
	likelihood=stackLikelihood,
  data = stack.dat, prior = stackPrior,
  sampler = stackSampler)

# Compare the sample distributions.
# Adjust for the multiple test scripts being run.
compareSampleDistributions(
  getSamples( bhlmSamples ), 
  getSamples( bugsSamples ), print.ks=F, 
  pvalue.cutoff = (0.01 / 100))
}
