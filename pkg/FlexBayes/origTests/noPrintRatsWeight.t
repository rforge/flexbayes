{cat("----- No extraneous info printed while fitting rats weight model --\n");T}
{
# Functions: bhlm, posteriorSamples
# Description: Make sure that no extraneous info is printed 
# while fitting a model to rat weights
  
  # save printed output to a file
  sink( "temp.txt" )
  
  # Load the data.  This is data for the weights of 
  # infant rats during the first 36 days of life.
  # The original study compares the growth rates of
  # the rats in a control group to those in a 
  # treatment group.  Here we estimate the growth 
  # rates in the control group only.
  data(ratsWeight)

  # Choose the hyperparameter values.  They
  # are taken to be equal to
  # those given in Gelfand et al. (1990).
  R <- diag( c( 100, 0.1 ) )
  rho <- 2

  # First, fit the model using bhlm.
  # Specify the priors.
  error.var.prior <- bayes.invChisq( df = 3, 
    sigma0.sq = 10 )
  random.var.prior <- bayes.invWishart( df = rho, 
    scale = solve( R ) )
  rats.prior <- bhlm.prior( error.var = error.var.prior,
    random.var = random.var.prior,
    level2.coef = "non-informative" )

  # Specify the control variates for the chain
  rats.sampler <- bhlm.sampler( nSamples = 1000,
    nThin = 10, init.point = "prior" )

  # Prepare the data to be fit using bhlm
  rats.bhlm <- data.frame( y = as.vector( ratsWeight$y ),
    x = rep( x = ratsWeight$x, each = 30 ),
    id = rep( (1:30), 5 ) )

  # Fit the model using bhlm
  bhlm.samples <- bhlm( data = rats.bhlm,
    random.formula = y ~ x, group.formula = ~id,
    level2.formula = ~ 1, 
    prior = rats.prior, sampler = rats.sampler )

  # Now fit a similar model using the posteriorSamples
  # function

  # Specify the BUGS model
  bugsModel <- function(){
    for ( i in 1:I ){
      for ( j in 1:J ){
        y[i,j] ~ dnorm(yMean[i,j], prec)
        yMean[i,j] <- beta[i,1] + beta[i,2] * x[j]
      }
      beta[i,1:2] ~ dmnorm (betaMean[], betaPrec[,])

    }
    betaMean[1:2] ~ dmnorm (eta[], Cinv[,])
    betaPrec[1:2,1:2] ~ dwish ( RtimesRho[,] , rho )
    betaCov[1:2,1:2] <- inverse( betaPrec[,] )
    RtimesRho[1,1] <- R[1,1] * rho
    RtimesRho[1,2] <- R[1,2] * rho
    RtimesRho[2,1] <- R[2,1] * rho
    RtimesRho[2,2] <- R[2,2] * rho
    prec ~ dgamma(precA, precB)
    sigma <- 1 / sqrt(prec)
    invisible()
  }

  # Prepare the data for the BUGS analysis
  rats.bugs <- ratsWeight
  rats.bugs$I <- dim(ratsWeight$y)[1]
  rats.bugs$J <- dim(ratsWeight$y)[2]

  # Specify the hyperparameter values
  rats.bugs$R <- R
  rats.bugs$rho <- rho
  rats.bugs$eta <- rep( 0, 2 )
  rats.bugs$Cinv <- diag( rep( 0.001, 2 ) )
  rats.bugs$precA <- 0.001
  rats.bugs$precB <- 0.001

  # Specify the initial values
  inits <- 
    list(prec = rats.bugs$precA / rats.bugs$precB,
      betaPrec = diag(2), betaMean = c(0,0))
  inits <- list(inits)

  # Select the parameters for which the posterior samples
  # are needed for inference.
  paramsToSave <- 
    c("betaCov", "betaMean", "sigma")

  bugs.samples <- posteriorSamples(model = bugsModel, 
    data = rats.bugs, inits = inits, 
    parametersToSave = paramsToSave, DIC = F, 
    engine = "WinBUGS")

  e <- effectiveSize(bugs.samples)
  a <- autocorr(bugs.samples)
  g <- geweke.diag(bugs.samples)
  h <- heidel.diag(bugs.samples)

  # The convergence diagnostics look good, so 
  # summarize the inferences
  s <- summary(bugs.samples)

  sink()
  ( file.lengths("temp.txt")[1] == 0 )
}
