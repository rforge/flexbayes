#########################################################################
##  
##  Validation of the bhlm function
##
##  
##

bhlmCompareBUGS <- function() {
	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

	# Fit the model using the bhlm function
	
	alpha.prior <- bayes.nonInformative()
	error.prior <- bayes.invChisq(3,1)
	betaCov.prior <- bayes.invChisq(3,1)
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  level2.coef = alpha.prior, 
	  random.var=betaCov.prior, common.error.var=2)
	
	orthoLike <- bhlm.likelihood(type="normal")
	orthoSampler <- bhlm.sampler( nBurnin=1000, 
	  nSamples=1000, nThin=50,
	  nChains=1, init.point = "prior")
	
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
	  level2.formula = ~Sex, group.formula = ~Subject, 
	  data = Orthodont, prior = orthoPrior, 
	  likelihood = orthoLike, sampler = orthoSampler)
	
	# Prepare the data set for the BUGS analysis
	Orthodont2 <- list(distance=Orthodont$distance,
		age=Orthodont$age)
	Nsubjects <- length(unique(Orthodont$Subject))
	# Create the sex vector to specify the sex for each
	# subject, rather than the sex associated with each
	# data point.
	# Also generate the numeric subject index.
	sex.bysubj <- rep(0, Nsubjects)
	Ndata <- length(Orthodont$distance)
	Orthodont2$Subject <- rep(0, Ndata)
	for (subj in (1:Nsubjects)){
		subjID <- unique(Orthodont$Subject)[subj]
		thisSubj <- which(subjID == Orthodont$Subject)
		Orthodont2$Subject[thisSubj] <- subj
		sex.bysubj[subj] <- unique(Orthodont$Sex[thisSubj] == "Female")
	}
	Orthodont2$Sex <- sex.bysubj
	Orthodont2$Nsubjects <- Nsubjects
	Orthodont2$Ndata <- Ndata
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], sigma2inv)
			mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
	  }
	  for (j in 1:Nsubjects){
	  	beta[j,1] ~ dnorm(mean.beta[1,j], prec.beta1)
	  	beta[j,2] ~ dnorm(mean.beta[2,j], prec.beta2)
	  	mean.beta[1,j] <- alpha00 + alpha01*Sex[j]
	  	mean.beta[2,j] <- alpha10 + alpha11*Sex[j]
	  }
	  prec.beta1 ~ dgamma(1.5, 1.5)
	  tau.beta1 <- 1/sqrt(prec.beta1)
	  prec.beta2 ~ dgamma(1.5, 1.5)
	  tau.beta2 <- 1/sqrt(prec.beta2)
	  sigma2inv ~ dgamma(1.5, 1.5)
	  sigma <- 1/sqrt(sigma2inv)
	  alpha00 ~ dnorm(0,1e-5)
		alpha01 ~ dnorm(0,1e-5)
		alpha10 ~ dnorm(0,1e-5)
		alpha11 ~ dnorm(0,1e-5)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("sigma", "alpha00", 
	  "alpha01", "alpha10", "alpha11", "tau.beta1", "tau.beta2",
	  "beta")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = Orthodont2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index = c(59, 1:4, 60:61, 5)
	bhlm.index = c((1:7), 36)
	compareSampleDistributions(getSamples(bugs.samples, bugs.index), 
	  getSamples(bhlm.samples, bhlm.index), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100)) 
}

bhlmCompareBUGS2 <- function() {
	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

	# Fit the model using the bhlm function
	
	alpha.prior <- bayes.normal(mean=zero, cov=identity)
	error.prior <- bayes.invChisq(3,1)
	betaCov.prior <- bayes.invChisq(3,1)
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  level2.coef = alpha.prior, 
	  random.var=betaCov.prior, common.error.var=2)
	
	orthoLike <- bhlm.likelihood(type="normal")
	
	orthoSampler <- bhlm.sampler( nBurnin=10000,
	  nSamples=1000, nThin=50,
	  nChains=1, init.point = "prior")
	
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
	  level2.formula = ~Sex, group.formula = ~Subject, 
	  data = Orthodont, prior = orthoPrior, 
	  likelihood = orthoLike, sampler = orthoSampler)
	
	# Prepare the data set for the BUGS analysis
	Orthodont2 <- list(distance=Orthodont$distance,
		age=Orthodont$age)
	Nsubjects <- length(unique(Orthodont$Subject))
	# Create the sex vector to specify the sex for each
	# subject, rather than the sex associated with each
	# data point.
	# Also generate the numeric subject index.
	sex.bysubj <- rep(0, Nsubjects)
	Ndata <- length(Orthodont$distance)
	Orthodont2$Subject <- rep(0, Ndata)
	for (subj in (1:Nsubjects)){
		subjID <- unique(Orthodont$Subject)[subj]
		thisSubj <- which(subjID == Orthodont$Subject)
		Orthodont2$Subject[thisSubj] <- subj
		sex.bysubj[subj] <- unique(Orthodont$Sex[thisSubj] == "Female")
	}
	Orthodont2$Sex <- sex.bysubj
	Orthodont2$Nsubjects <- Nsubjects
	Orthodont2$Ndata <- Ndata
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], sigma2inv)
			mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
	  }
	  for (j in 1:Nsubjects){
	  	beta[j,1] ~ dnorm(mean.beta[1,j], prec.beta1)
	  	beta[j,2] ~ dnorm(mean.beta[2,j], prec.beta2)
	  	mean.beta[1,j] <- alpha00 + alpha01*Sex[j]
	  	mean.beta[2,j] <- alpha10 + alpha11*Sex[j]
	  }
	  prec.beta1 ~ dgamma(1.5, 1.5)
	  tau.beta1 <- 1/sqrt(prec.beta1)
	  prec.beta2 ~ dgamma(1.5, 1.5)
	  tau.beta2 <- 1/sqrt(prec.beta2)
	  sigma2inv ~ dgamma(1.5, 1.5)
	  sigma <- 1/sqrt(sigma2inv)
	  alpha00 ~ dnorm(0,1)
		alpha01 ~ dnorm(0,1)
		alpha10 ~ dnorm(0,1)
		alpha11 ~ dnorm(0,1)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("sigma", "alpha00", 
	  "alpha01", "alpha10", "alpha11", "tau.beta1", "tau.beta2", 
	  "beta")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = Orthodont2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index = c(59, 1:4, 60:61, 5)
	bhlm.index = c((1:7), 36)
	compareSampleDistributions(getSamples(bugs.samples, bugs.index), 
	  getSamples(bhlm.samples, bhlm.index), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100)) 
}

bhlmCompareBUGS3 <- function() {

	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

	# Fit the model using the bhlm function
	
	ortho <- as.data.frame( 
	  Orthodont[ Orthodont$Sex == "Female", ] )
	
	alpha.prior <- bayes.normal(mean=zero, cov=identity)
	error.prior <- bayes.invChisq(3,1)
	randomvar.prior <- bayes.invChisq(3,1)
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  level2.coef = alpha.prior, 
	  random.var=randomvar.prior, common.error.var=2)
	
	orthoLike <- bhlm.likelihood(type="normal")
	
	orthoSampler <- bhlm.sampler( nBurnin=10000, 
	  nSamples=1000, nThin=50,
	  nChains=1, init.point = "prior")
	
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
	  group.formula = ~Subject, level2.formula = ~ 1,
	  data = ortho, prior = orthoPrior, 
	  likelihood = orthoLike, sampler = orthoSampler)
	
	# Prepare the data set for the BUGS analysis
	ortho2 <- list(distance=ortho$distance,
		age=ortho$age)
	Nsubjects <- length(unique(ortho$Subject))
	Ndata <- length(ortho$distance)
	# generate the numeric subject index.
	for (subj in (1:Nsubjects)){
		subjID <- unique(ortho$Subject)[subj]
		thisSubj <- which(subjID == ortho$Subject)
		ortho2$Subject[thisSubj] <- subj
	}
	ortho2$Nsubjects <- Nsubjects
	ortho2$Ndata <- Ndata
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], sigma2inv)
			mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
	  }
	  for (j in 1:Nsubjects){
	  	beta[j,1] ~ dnorm(mean.beta[1,j], prec.beta1)
	  	beta[j,2] ~ dnorm(mean.beta[2,j], prec.beta2)
	  	mean.beta[1,j] <- alpha0
	  	mean.beta[2,j] <- alpha1
	  }
	  prec.beta1 ~ dgamma(1.5, 1.5)
	  tau.beta1 <- 1/sqrt(prec.beta1)
	  prec.beta2 ~ dgamma(1.5, 1.5)
	  tau.beta2 <- 1/sqrt(prec.beta2)
	  sigma2inv ~ dgamma(1.5, 1.5)
	  sigma <- 1/sqrt(sigma2inv)
	  alpha0 ~ dnorm(0,1)
		alpha1 ~ dnorm(0,1)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("sigma", "alpha0", 
		"alpha1", "tau.beta1", "tau.beta2", "beta")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = ortho2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index <- c(25, 1:2, 26:27, 3)
	bhlm.index <- c((1:5), 12)
	compareSampleDistributions(getSamples(bugs.samples, bugs.index), 
	  getSamples(bhlm.samples, bhlm.index), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100))
}

bhlmCompareBUGS4 <- function() {
	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

	# Fit the model using the bhlm function
	
	alpha.prior <- bayes.normal(rep(0,4), diag(rep(1,4)))
	error.prior <- bayes.nonInfoPower(-1)
	betaCov.prior <- bayes.invChisq(3,1)
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  level2.coef = alpha.prior, 
	  random.var=betaCov.prior, common.error.var=2)
	
	orthoLike <- bhlm.likelihood(type="normal")
	
	orthoSampler <- bhlm.sampler( nBurnin=20000, 
	  nSamples=1000, nThin=50,
	  nChains=1, init.point = "prior")
	
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
	  level2.formula = ~Sex, group.formula = ~Subject, 
	  data = Orthodont, prior = orthoPrior, 
	  likelihood = orthoLike, sampler = orthoSampler)
	
	# Prepare the data set for the BUGS analysis
	Orthodont2 <- list(distance=Orthodont$distance,
		age=Orthodont$age)
	Nsubjects <- length(unique(Orthodont$Subject))
	# Create the sex vector to specify the sex for each
	# subject, rather than the sex associated with each
	# data point.
	# Also generate the numeric subject index.
	sex.bysubj <- rep(0, Nsubjects)
	Ndata <- length(Orthodont$distance)
	Orthodont2$Subject <- rep(0, Ndata)
	for (subj in (1:Nsubjects)){
		subjID <- unique(Orthodont$Subject)[subj]
		thisSubj <- which(subjID == Orthodont$Subject)
		Orthodont2$Subject[thisSubj] <- subj
		sex.bysubj[subj] <- unique(Orthodont$Sex[thisSubj] == "Female")
	}
	Orthodont2$Sex <- sex.bysubj
	Orthodont2$Nsubjects <- Nsubjects
	Orthodont2$Ndata <- Ndata
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], sigma2inv)
			mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
	  }
	  for (j in 1:Nsubjects){
	  	beta[j,1] ~ dnorm(mean.beta[1,j], prec.beta1)
	  	beta[j,2] ~ dnorm(mean.beta[2,j], prec.beta2)
	  	mean.beta[1,j] <- alpha00 + alpha01*Sex[j]
	  	mean.beta[2,j] <- alpha10 + alpha11*Sex[j]
	  }
	  prec.beta1 ~ dgamma(1.5, 1.5)
	  tau.beta1 <- 1/sqrt(prec.beta1)
	  prec.beta2 ~ dgamma(1.5, 1.5)
	  tau.beta2 <- 1/sqrt(prec.beta2)
	  sigma2inv ~ dgamma(1e-2, 1e-6)
	  sigma <- 1/sqrt(sigma2inv)
	  alpha00 ~ dnorm(0,1)
		alpha01 ~ dnorm(0,1)
		alpha10 ~ dnorm(0,1)
		alpha11 ~ dnorm(0,1)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("beta", "alpha00", "alpha01", 
		"alpha10", "alpha11", "sigma", "tau.beta1", "tau.beta2")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = Orthodont2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# compare the sample distributions for the alpha parameters
	# and tau parameters; do not do so for the other parameters. The 
	# posterior inferences are close but not exact,
	# due to the prior approximation here
	bugs.index <- c( 59, (1:4), 60:61, (5:58) )
	bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples, (2:7)), 
	  getSamples(bhlm.samples, (2:7)), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100))
}

bhlmCompareBUGS5 <- function() {

	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

	# Fit the model using the bhlm function
	
	alpha.prior <- bayes.nonInformative()
	error.prior <- bayes.invChisq(3,1)
	betaCov.prior <- bayes.invChisq(3,1)
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  level2.coef = alpha.prior, 
	  random.var=betaCov.prior, common.error.var=2)
	
	orthoLike <- bhlm.likelihood(type="normal")
	orthoSampler <- bhlm.sampler( nBurnin=10000, 
	  nSamples=1000, nThin=50,
	  nChains=1, init.point = "prior")
	
	bhlm.samples <- bhlm(random.formula = distance~1,
	  level2.formula = ~Sex, group.formula = ~Subject, 
	  data = Orthodont, prior = orthoPrior, 
	  likelihood = orthoLike, sampler = orthoSampler)
	
	# Prepare the data set for the BUGS analysis
	Orthodont2 <- list(distance=Orthodont$distance)
	Nsubjects <- length(unique(Orthodont$Subject))
	# Create the sex vector to specify the sex for each
	# subject, rather than the sex associated with each
	# data point.
	# Also generate the numeric subject index.
	sex.bysubj <- rep(0, Nsubjects)
	Ndata <- length(Orthodont$distance)
	Orthodont2$Subject <- rep(0, Ndata)
	for (subj in (1:Nsubjects)){
		subjID <- unique(Orthodont$Subject)[subj]
		thisSubj <- which(subjID == Orthodont$Subject)
		Orthodont2$Subject[thisSubj] <- subj
		sex.bysubj[subj] <- unique(Orthodont$Sex[thisSubj] == "Female")
	}
	Orthodont2$Sex <- sex.bysubj
	Orthodont2$Nsubjects <- Nsubjects
	Orthodont2$Ndata <- Ndata
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			distance[i] ~ dnorm(mean[i], sigma2inv)
			mean[i] <- beta[Subject[i],1]
	  }
	  for (j in 1:Nsubjects){
	  	beta[j,1] ~ dnorm(mean.beta[1,j], prec.beta1)
	  	mean.beta[1,j] <- alpha00 + alpha01*Sex[j]
	  }
	  prec.beta1 ~ dgamma(1.5, 1.5)
	  tau.beta1 <- 1/sqrt(prec.beta1)
	  sigma2inv ~ dgamma(1.5, 1.5)
	  sigma <- 1/sqrt(sigma2inv)
	  alpha00 ~ dnorm(0,1e-5)
		alpha01 ~ dnorm(0,1e-5)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("sigma", "alpha00", 
	  "alpha01", "tau.beta1", "beta")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = Orthodont2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# Compare the sample distributions.
	bugs.index = c(30, 1:2, 31, 3:4, 19)
	bhlm.index = c((1:4), 19, 7, 24)
	# Adjust for the multiple test scripts being run.
	compareSampleDistributions(getSamples(bugs.samples, bugs.index), 
	  getSamples(bhlm.samples, bhlm.index), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100)) 
}

bhlmCompareBUGS6 <- function() {
	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

	# Fit the model using the bhlm function
	
	error.prior <- bayes.invChisq(3,1)
  random.var.prior = bayes.invWishart( df = 3, scale = diag( c(2,2) ) )
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  random.var=random.var.prior)
	
	orthoLike <- bhlm.likelihood(type="normal")
	
	orthoSampler <- bhlm.sampler( nBurnin=10000, 
	  nSamples=1000, nThin=100,
	  nChains=1, init.point = "prior")
	
	orthoDat <- Orthodont[1:4,]
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
	  data = orthoDat, prior = orthoPrior, 
	  likelihood = orthoLike, sampler = orthoSampler, 
	  debug = F)
	
	# Prepare the data set for the BUGS analysis
	orthoDat2 <- list(distance=orthoDat$distance,
		age=orthoDat$age)
	Ndata <- length(orthoDat$distance)
  orthoDat2$Ndata <- Ndata
	orthoDat2$tau2inv.prec = diag( c(0.5,0.5) ) 
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], sigma2inv)
			mean[i] <- beta[1] + beta[2]*age.centered[i]
	  }
	  beta[1:2] ~ dmnorm( mean.beta[], tau2inv[,])
	  mean.beta[1] <- 0
	  mean.beta[2] <- 0
	  
	 	tau2inv[1:2,1:2] ~ dwish( tau2inv.prec[,], 3 )
	 	tau2[1:2,1:2] <- inverse( tau2inv[,] )
	  sigma2inv ~ dgamma(1.5, 1.5)
	  sigma <- 1/sqrt(sigma2inv)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("beta", "sigma", "tau2")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = orthoDat2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# compare the sample distributions 
	bugs.index <- c( 3:7, 1:2 )
	bugs.samples <- bugs.samples[,bugs.index]
	compareSampleDistributions(getSamples(bugs.samples), 
	  getSamples(bhlm.samples), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100))
}

bhlmCompareBUGS7 <- function() {

	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

  ## Like bhlmCompareBUGS2, with different prior hyperparameters.
  
	# Fit the model using the bhlm function
	
	alpha.prior <- bayes.normal(mean=rep(2,4), cov=rep(1/5.5, 4))
	error.prior <- bayes.invChisq(6,0.5)
	betaCov.prior <- bayes.invChisq(5,2)
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  level2.coef = alpha.prior, 
	  random.var=betaCov.prior, common.error.var=2)
	
	orthoLike <- bhlm.likelihood(type="normal")
	
	orthoSampler <- bhlm.sampler( nBurnin=10000,
	  nSamples=1000, nThin=50,
	  nChains=1, init.point = "prior")
	
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
	  level2.formula = ~Sex, group.formula = ~Subject, 
	  data = Orthodont, prior = orthoPrior, 
	  likelihood = orthoLike, sampler = orthoSampler)
	
	# Prepare the data set for the BUGS analysis
	Orthodont2 <- list(distance=Orthodont$distance,
		age=Orthodont$age)
	Nsubjects <- length(unique(Orthodont$Subject))
	# Create the sex vector to specify the sex for each
	# subject, rather than the sex associated with each
	# data point.
	# Also generate the numeric subject index.
	sex.bysubj <- rep(0, Nsubjects)
	Ndata <- length(Orthodont$distance)
	Orthodont2$Subject <- rep(0, Ndata)
	for (subj in (1:Nsubjects)){
		subjID <- unique(Orthodont$Subject)[subj]
		thisSubj <- which(subjID == Orthodont$Subject)
		Orthodont2$Subject[thisSubj] <- subj
		sex.bysubj[subj] <- unique(Orthodont$Sex[thisSubj] == "Female")
	}
	Orthodont2$Sex <- sex.bysubj
	Orthodont2$Nsubjects <- Nsubjects
	Orthodont2$Ndata <- Ndata
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], sigma2inv)
			mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
	  }
	  for (j in 1:Nsubjects){
	  	beta[j,1] ~ dnorm(mean.beta[1,j], prec.beta1)
	  	beta[j,2] ~ dnorm(mean.beta[2,j], prec.beta2)
	  	mean.beta[1,j] <- alpha00 + alpha01*Sex[j]
	  	mean.beta[2,j] <- alpha10 + alpha11*Sex[j]
	  }
	  prec.beta1 ~ dgamma(2.5, 5)
	  tau.beta1 <- 1/sqrt(prec.beta1)
	  prec.beta2 ~ dgamma(2.5, 5)
	  tau.beta2 <- 1/sqrt(prec.beta2)
	  sigma2inv ~ dgamma(3, 1.5)
	  sigma <- 1/sqrt(sigma2inv)
	  alpha00 ~ dnorm(2,5.5)
		alpha01 ~ dnorm(2,5.5)
		alpha10 ~ dnorm(2,5.5)
		alpha11 ~ dnorm(2,5.5)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("sigma", "alpha00", 
	  "alpha01", "alpha10", "alpha11", "tau.beta1", "tau.beta2", 
	  "beta")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = Orthodont2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index = c(59, 1:4, 60:61, 5)
	bhlm.index = c((1:7), 36)
	compareSampleDistributions(getSamples(bugs.samples, bugs.index), 
	  getSamples(bhlm.samples, bhlm.index), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100)) 
}

bhlmCompareBUGS8 <- function() {

	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

  ## "Non-informative" power law prior for the random effect variances
  ## The results only match for very few parameters, since it is not 
  ## possible to match the priors exactly.
  
	# Fit the model using the bhlm function
	
	alpha.prior <- bayes.normal(mean=rep(2,4), cov=rep(1/5.5, 4))
	error.prior <- bayes.invChisq(6,0.5)
	betaCov.prior <- bayes.nonInfoPower(-1)
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  level2.coef = alpha.prior, 
	  random.var=betaCov.prior, common.error.var=2)
	
	orthoLike <- bhlm.likelihood(type="normal")
	
	orthoSampler <- bhlm.sampler( nBurnin=10000,
	  nSamples=1000, nThin=50,
	  nChains=1, init.point = "prior")
	
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
	  level2.formula = ~Sex, group.formula = ~Subject, 
	  data = Orthodont, prior = orthoPrior, 
	  likelihood = orthoLike, sampler = orthoSampler)
	
	# Prepare the data set for the BUGS analysis
	Orthodont2 <- list(distance=Orthodont$distance,
		age=Orthodont$age)
	Nsubjects <- length(unique(Orthodont$Subject))
	# Create the sex vector to specify the sex for each
	# subject, rather than the sex associated with each
	# data point.
	# Also generate the numeric subject index.
	sex.bysubj <- rep(0, Nsubjects)
	Ndata <- length(Orthodont$distance)
	Orthodont2$Subject <- rep(0, Ndata)
	for (subj in (1:Nsubjects)){
		subjID <- unique(Orthodont$Subject)[subj]
		thisSubj <- which(subjID == Orthodont$Subject)
		Orthodont2$Subject[thisSubj] <- subj
		sex.bysubj[subj] <- unique(Orthodont$Sex[thisSubj] == "Female")
	}
	Orthodont2$Sex <- sex.bysubj
	Orthodont2$Nsubjects <- Nsubjects
	Orthodont2$Ndata <- Ndata
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], sigma2inv)
			mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
	  }
	  for (j in 1:Nsubjects){
	  	beta[j,1] ~ dnorm(mean.beta[1,j], prec.beta1)
	  	beta[j,2] ~ dnorm(mean.beta[2,j], prec.beta2)
	  	mean.beta[1,j] <- alpha00 + alpha01*Sex[j]
	  	mean.beta[2,j] <- alpha10 + alpha11*Sex[j]
	  }
	  prec.beta1 ~ dgamma(1e-3, 1e-6)
	  tau.beta1 <- 1/sqrt(prec.beta1)
	  prec.beta2 ~ dgamma(1e-3, 1e-6)
	  tau.beta2 <- 1/sqrt(prec.beta2)
	  sigma2inv ~ dgamma(3, 1.5)
	  sigma <- 1/sqrt(sigma2inv)
	  alpha00 ~ dnorm(2,5.5)
		alpha01 ~ dnorm(2,5.5)
		alpha10 ~ dnorm(2,5.5)
		alpha11 ~ dnorm(2,5.5)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("sigma", "alpha00", 
	  "alpha01", "alpha10", "alpha11", "tau.beta1", "tau.beta2", 
	  "beta")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = Orthodont2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index = c(59, 1:4, 60:61, 5)
	bhlm.index = c((1:7), 36)
	compareSampleDistributions(getSamples(bugs.samples, bugs.index), 
	  getSamples(bhlm.samples, bhlm.index), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100)) 
}

bhlmCompareBUGS9 <- function() {
	
	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

	# Fit the model using the bhlm function
	
	error.prior <- bayes.invChisq(3,1)
  random.var.prior = bayes.invWishart( df = 3, scale = diag( c(2,2) ) )
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  random.var=random.var.prior)
	
	orthoSampler <- bhlm.sampler( nBurnin=10000, 
	  nSamples=1000, nThin=100,
	  nChains=1, init.point = "prior")
	
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
    group.formula = ~Subject, 
	  data = Orthodont, prior = orthoPrior, 
	  sampler = orthoSampler, 
	  debug = F)
	
	# Prepare the data set for the BUGS analysis
	Orthodont2 <- list(distance=Orthodont$distance,
		age=Orthodont$age)
	Nsubjects <- length(unique(Orthodont$Subject))
	Ndata <- length(Orthodont$distance)
	Orthodont2$Subject <- rep(0, Ndata)
	for (subj in (1:Nsubjects)){
		subjID <- unique(Orthodont$Subject)[subj]
		thisSubj <- which(subjID == Orthodont$Subject)
		Orthodont2$Subject[thisSubj] <- subj
	}
	Orthodont2$Nsubjects <- Nsubjects
	Orthodont2$Ndata <- Ndata
	Orthodont2$tau2inv.prec = diag( c(0.5,0.5) ) 
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], sigma2inv)
			mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
	  }
	  for (j in 1:Nsubjects){
	    beta[j,1:2] ~ dmnorm( mean.beta[], tau2inv[,])	
	  }
	  mean.beta[1] <- 0
	  mean.beta[2] <- 0
	 	tau2inv[1:2,1:2] ~ dwish( tau2inv.prec[,], 3 )
	 	tau2[1:2,1:2] <- inverse( tau2inv[,] )
	  sigma2inv ~ dgamma(1.5, 1.5)
	  sigma <- 1/sqrt(sigma2inv)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("beta", "sigma", "tau2")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = Orthodont2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# compare the sample distributions 
	bugs.index <- c( 55:59, 1:2 )
	bhlm.index <- c( 1:5, 34:35 )
	compareSampleDistributions(
	  getSamples(bugs.samples[,bugs.index]), 
	  getSamples(bhlm.samples[,bhlm.index]), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100))
}

bhlmCompareBUGS10 <- function() {
	
	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

	# Fit the model using the bhlm function
	
	alpha.prior <- bayes.normal(rep(0,4), diag(rep(1,4)))
	error.prior <- bayes.invChisq(3,1)
  random.var.prior = bayes.invWishart( df = 3, scale = diag( c(2,2) ) )
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  level2.coef = alpha.prior, 
	  random.var=random.var.prior)
	
	orthoSampler <- bhlm.sampler( nBurnin=10000, 
	  nSamples=1000, nThin=100,
	  nChains=1, init.point = "prior")
	
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
	  level2.formula = ~Sex, group.formula = ~Subject, 
	  data = Orthodont, prior = orthoPrior, 
	  sampler = orthoSampler, 
	  debug = F)
	
	# Prepare the data set for the BUGS analysis
	Orthodont2 <- list(distance=Orthodont$distance,
		age=Orthodont$age)
	Nsubjects <- length(unique(Orthodont$Subject))
	# Create the sex vector to specify the sex for each
	# subject, rather than the sex associated with each
	# data point.
	# Also generate the numeric subject index.
	sex.bysubj <- rep(0, Nsubjects)
	Ndata <- length(Orthodont$distance)
	Orthodont2$Subject <- rep(0, Ndata)
	for (subj in (1:Nsubjects)){
		subjID <- unique(Orthodont$Subject)[subj]
		thisSubj <- which(subjID == Orthodont$Subject)
		Orthodont2$Subject[thisSubj] <- subj
		sex.bysubj[subj] <- unique(Orthodont$Sex[thisSubj] == "Female")
	}
	Orthodont2$Sex <- sex.bysubj
	Orthodont2$Nsubjects <- Nsubjects
	Orthodont2$Ndata <- Ndata
	Orthodont2$tau2inv.prec = diag( c(0.5,0.5) ) 
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], sigma2inv)
			mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
	  }
	  for (j in 1:Nsubjects){
	    beta[j,1:2] ~ dmnorm( mean.beta[j,], tau2inv[,])
	  	mean.beta[j,1] <- alpha00 + alpha01*Sex[j]
	  	mean.beta[j,2] <- alpha10 + alpha11*Sex[j]
	  }
	 	tau2inv[1:2,1:2] ~ dwish( tau2inv.prec[,], 3 )
	 	tau2[1:2,1:2] <- inverse( tau2inv[,] )
	  sigma2inv ~ dgamma(1.5, 1.5)
	  sigma <- 1/sqrt(sigma2inv)
	  alpha00 ~ dnorm(0,1)
		alpha01 ~ dnorm(0,1)
		alpha10 ~ dnorm(0,1)
		alpha11 ~ dnorm(0,1)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("beta", "alpha00", "alpha01", 
		"alpha10", "alpha11", "sigma", "tau2")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = Orthodont2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# compare the sample distributions 
	bugs.index <- c( 59, (1:4), 60:63, (5:6) )
  bhlm.index <- c( 1:9, 38:39 )
	compareSampleDistributions(
	  getSamples(bugs.samples[,bugs.index]), 
	  getSamples(bhlm.samples[,bhlm.index]), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100))
}

bhlmCompareBUGS11 <- function() {

	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

  # test manual specification of initial values
  
	# Fit the model using the bhlm function
	
	alpha.prior <- bayes.nonInformative()
	error.prior <- bayes.invChisq(3,1)
	betaCov.prior <- bayes.invChisq(3,1)
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  level2.coef = alpha.prior, 
	  random.var=betaCov.prior)
	
	orthoLike <- bhlm.likelihood(type="normal")
	orthoSampler <- bhlm.sampler( nBurnin=1000, 
	  nSamples=1000, nThin=50,
	  nChains=1, init.point = "user's choice",
	  params.init = list( error.var = 1.0, 
	  random.coef = matrix( 0.0, nrow = 2, ncol = 27 ),
	  level2.coef = rep(0, 4), random.var = c(1,1),
	  fixed.coef = 0 ) )
	
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
	  level2.formula = ~Sex, group.formula = ~Subject, 
	  data = Orthodont, prior = orthoPrior, 
	  likelihood = orthoLike, sampler = orthoSampler)
	
	# Prepare the data set for the BUGS analysis
	Orthodont2 <- list(distance=Orthodont$distance,
		age=Orthodont$age)
	Nsubjects <- length(unique(Orthodont$Subject))
	# Create the sex vector to specify the sex for each
	# subject, rather than the sex associated with each
	# data point.
	# Also generate the numeric subject index.
	sex.bysubj <- rep(0, Nsubjects)
	Ndata <- length(Orthodont$distance)
	Orthodont2$Subject <- rep(0, Ndata)
	for (subj in (1:Nsubjects)){
		subjID <- unique(Orthodont$Subject)[subj]
		thisSubj <- which(subjID == Orthodont$Subject)
		Orthodont2$Subject[thisSubj] <- subj
		sex.bysubj[subj] <- unique(Orthodont$Sex[thisSubj] == "Female")
	}
	Orthodont2$Sex <- sex.bysubj
	Orthodont2$Nsubjects <- Nsubjects
	Orthodont2$Ndata <- Ndata
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], sigma2inv)
			mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
	  }
	  for (j in 1:Nsubjects){
	  	beta[j,1] ~ dnorm(mean.beta[1,j], prec.beta1)
	  	beta[j,2] ~ dnorm(mean.beta[2,j], prec.beta2)
	  	mean.beta[1,j] <- alpha00 + alpha01*Sex[j]
	  	mean.beta[2,j] <- alpha10 + alpha11*Sex[j]
	  }
	  prec.beta1 ~ dgamma(1.5, 1.5)
	  tau.beta1 <- 1/sqrt(prec.beta1)
	  prec.beta2 ~ dgamma(1.5, 1.5)
	  tau.beta2 <- 1/sqrt(prec.beta2)
	  sigma2inv ~ dgamma(1.5, 1.5)
	  sigma <- 1/sqrt(sigma2inv)
	  alpha00 ~ dnorm(0,1e-5)
		alpha01 ~ dnorm(0,1e-5)
		alpha10 ~ dnorm(0,1e-5)
		alpha11 ~ dnorm(0,1e-5)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("sigma", "alpha00", 
	  "alpha01", "alpha10", "alpha11", "tau.beta1", "tau.beta2",
	  "beta")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = Orthodont2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index = c(59, 1:4, 60:61, 5)
	bhlm.index = c((1:7), 36)
	compareSampleDistributions(getSamples(bugs.samples, bugs.index), 
	  getSamples(bhlm.samples, bhlm.index), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100)) 
}

bhlmCompareBUGSfixedErrorVar <- function() {
	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

	# Fit the model using the bhlm function
	
	alpha.prior <- bayes.nonInformative()
	error.prior <- bayes.massPoint(1.5)
	betaCov.prior <- bayes.invChisq(3,1)
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  level2.coef = alpha.prior, 
	  random.var=betaCov.prior, common.error.var=2)
	
	orthoLike <- bhlm.likelihood(type="normal")
	
	orthoSampler <- bhlm.sampler( nBurnin=1000, 
	  nSamples=1000, nThin=50,
	  nChains=1, init.point = "prior")
	
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
	  level2.formula = ~Sex, group.formula = ~Subject, 
	  data = Orthodont, prior = orthoPrior, 
	  likelihood = orthoLike, sampler = orthoSampler)
	
	# Prepare the data set for the BUGS analysis
	Orthodont2 <- list(distance=Orthodont$distance,
		age=Orthodont$age)
	Nsubjects <- length(unique(Orthodont$Subject))
	# Create the sex vector to specify the sex for each
	# subject, rather than the sex associated with each
	# data point.
	# Also generate the numeric subject index.
	sex.bysubj <- rep(0, Nsubjects)
	Ndata <- length(Orthodont$distance)
	Orthodont2$Subject <- rep(0, Ndata)
	for (subj in (1:Nsubjects)){
		subjID <- unique(Orthodont$Subject)[subj]
		thisSubj <- which(subjID == Orthodont$Subject)
		Orthodont2$Subject[thisSubj] <- subj
		sex.bysubj[subj] <- unique(Orthodont$Sex[thisSubj] == "Female")
	}
	Orthodont2$Sex <- sex.bysubj
	Orthodont2$Nsubjects <- Nsubjects
	Orthodont2$Ndata <- Ndata
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], 0.6667)
			mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
	  }
	  for (j in 1:Nsubjects){
	  	beta[j,1] ~ dnorm(mean.beta[1,j], prec.beta1)
	  	beta[j,2] ~ dnorm(mean.beta[2,j], prec.beta2)
	  	mean.beta[1,j] <- alpha00 + alpha01*Sex[j]
	  	mean.beta[2,j] <- alpha10 + alpha11*Sex[j]
	  }
	  prec.beta1 ~ dgamma(1.5, 1.5)
	  tau.beta1 <- 1/sqrt(prec.beta1)
	  prec.beta2 ~ dgamma(1.5, 1.5)
	  tau.beta2 <- 1/sqrt(prec.beta2)
	  alpha00 ~ dnorm(0,1e-5)
		alpha01 ~ dnorm(0,1e-5)
		alpha10 ~ dnorm(0,1e-5)
		alpha11 ~ dnorm(0,1e-5)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("alpha00", "alpha01", 
		"alpha10", "alpha11", "tau.beta1", "tau.beta2",
		"beta")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = Orthodont2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index <- c(1:4, 59:60, 5)
	bhlm.index <- c((1:6), 35)
	compareSampleDistributions(getSamples(bugs.samples, bugs.index), 
	  getSamples(bhlm.samples, bhlm.index), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100))
}

bhlmCompareBUGSfixedErrorVec <- function() {
	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

	# Fit the model using the bhlm function
	
	alpha.prior <- bayes.nonInformative()
	# fix the error variance for each group (subject)
	error.prior <- vector( "list", 27 )
	for( i in (1:27) ){
	  error.prior[[i]] <- bayes.massPoint(1.5)
	}
	betaCov.prior <- bayes.invChisq(3,1)
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  level2.coef = alpha.prior, 
	  random.var=betaCov.prior, common.error.var=0)
	
	orthoLike <- bhlm.likelihood(type="normal")
	
	orthoSampler <- bhlm.sampler( nBurnin=1000, 
	  nSamples=1000, nThin=50,
	  nChains=1, init.point = "prior")
	
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
	  level2.formula = ~Sex, group.formula = ~Subject, 
	  data = Orthodont, prior = orthoPrior, 
	  likelihood = orthoLike, sampler = orthoSampler)
	
	# Prepare the data set for the BUGS analysis
	Orthodont2 <- list(distance=Orthodont$distance,
		age=Orthodont$age)
	Nsubjects <- length(unique(Orthodont$Subject))
	# Create the sex vector to specify the sex for each
	# subject, rather than the sex associated with each
	# data point.
	# Also generate the numeric subject index.
	sex.bysubj <- rep(0, Nsubjects)
	Ndata <- length(Orthodont$distance)
	Orthodont2$Subject <- rep(0, Ndata)
	for (subj in (1:Nsubjects)){
		subjID <- unique(Orthodont$Subject)[subj]
		thisSubj <- which(subjID == Orthodont$Subject)
		Orthodont2$Subject[thisSubj] <- subj
		sex.bysubj[subj] <- unique(Orthodont$Sex[thisSubj] == "Female")
	}
	Orthodont2$Sex <- sex.bysubj
	Orthodont2$Nsubjects <- Nsubjects
	Orthodont2$Ndata <- Ndata
	
	## specify the model in BUGS format	
	 bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], 0.666667)
			mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
	  }
	  for (j in 1:Nsubjects){
	  	beta[j,1] ~ dnorm(mean.beta[1,j], prec.beta1)
	  	beta[j,2] ~ dnorm(mean.beta[2,j], prec.beta2)
	  	mean.beta[1,j] <- alpha00 + alpha01*Sex[j]
	  	mean.beta[2,j] <- alpha10 + alpha11*Sex[j]
	  }
	  prec.beta1 ~ dgamma(1.5, 1.5)
	  tau.beta1 <- 1/sqrt(prec.beta1)
	  prec.beta2 ~ dgamma(1.5, 1.5)
	  tau.beta2 <- 1/sqrt(prec.beta2)
	  alpha00 ~ dnorm(0,1e-6)
		alpha01 ~ dnorm(0,1e-6)
		alpha10 ~ dnorm(0,1e-6)
		alpha11 ~ dnorm(0,1e-6)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("alpha00", "alpha01", 
		"alpha10", "alpha11", "tau.beta1", "tau.beta2",
		"beta")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = Orthodont2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# Compare the sample distributions.
	# Adjust for the multiple test scripts being run.
	bugs.index <- c(1:4, 59:60, 5)
	bhlm.index <- c((1:6), 35)
	compareSampleDistributions(getSamples(bugs.samples, bugs.index), 
	  getSamples(bhlm.samples, bhlm.index), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100))
}

bhlmCompareBUGSnonInfo <- function() {
	# Set default contrasts so test comparisons work:
	oldOpt <- options(contrasts=c("contr.treatment", "contr.poly"))
        on.exit(options(oldOpt))

	# Fit the model using the bhlm function
	
	alpha.prior <- bayes.nonInformative()
	error.prior <- bayes.nonInfoPower(-1)
	betaCov.prior <- bayes.invChisq(3,1)
	
	orthoPrior <- bhlm.prior(error.var=error.prior, 
	  level2.coef = alpha.prior, 
	  random.var=betaCov.prior, common.error.var=2)
	
	orthoLike <- bhlm.likelihood(type="normal")
	
	orthoSampler <- bhlm.sampler( nBurnin=1000, 
	  nSamples=1000, nThin=50,
	  nChains=1, init.point = "prior")
	
	bhlm.samples <- bhlm(random.formula = distance~I(age-11),
	  level2.formula = ~Sex, group.formula = ~Subject, 
	  data = Orthodont, prior = orthoPrior, 
	  likelihood = orthoLike, sampler = orthoSampler)
	
	# Prepare the data set for the BUGS analysis
	Orthodont2 <- list(distance=Orthodont$distance,
		age=Orthodont$age)
	Nsubjects <- length(unique(Orthodont$Subject))
	# Create the sex vector to specify the sex for each
	# subject, rather than the sex associated with each
	# data point.
	# Also generate the numeric subject index.
	sex.bysubj <- rep(0, Nsubjects)
	Ndata <- length(Orthodont$distance)
	Orthodont2$Subject <- rep(0, Ndata)
	for (subj in (1:Nsubjects)){
		subjID <- unique(Orthodont$Subject)[subj]
		thisSubj <- which(subjID == Orthodont$Subject)
		Orthodont2$Subject[thisSubj] <- subj
		sex.bysubj[subj] <- unique(Orthodont$Sex[thisSubj] == "Female")
	}
	Orthodont2$Sex <- sex.bysubj
	Orthodont2$Nsubjects <- Nsubjects
	Orthodont2$Ndata <- Ndata
	
	## specify the model in BUGS format	
	bugsModel <- function (){
	  for( i in 1 : Ndata ) {
			age.centered[i] <- age[i] - 11
			distance[i] ~ dnorm(mean[i], sigma2inv)
			mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
	  }
	  for (j in 1:Nsubjects){
	  	beta[j,1] ~ dnorm(mean.beta[1,j], prec.beta1)
	  	beta[j,2] ~ dnorm(mean.beta[2,j], prec.beta2)
	  	mean.beta[1,j] <- alpha00 + alpha01*Sex[j]
	  	mean.beta[2,j] <- alpha10 + alpha11*Sex[j]
	  }
	  prec.beta1 ~ dgamma(1.5, 1.5)
	  tau.beta1 <- 1/sqrt(prec.beta1)
	  prec.beta2 ~ dgamma(1.5, 1.5)
	  tau.beta2 <- 1/sqrt(prec.beta2)
	  sigma2inv ~ dgamma(1e-2, 1e-6)
	  sigma <- 1/sqrt(sigma2inv)
	  alpha00 ~ dnorm(0,1e-5)
		alpha01 ~ dnorm(0,1e-5)
		alpha10 ~ dnorm(0,1e-5)
		alpha11 ~ dnorm(0,1e-5)
		
	  invisible()		
	}
	
	# Specify the parameters for which we wish to save 
	# the posterior samples
	parameters.to.save <- c("sigma", "alpha00", "alpha01", 
		"alpha10", "alpha11", "tau.beta1", "tau.beta2",
		"beta")
	
	## obtain the posterior samples
	bugs.samples <- posteriorSamples (
	  data = Orthodont2, model = bugsModel, 
	  nChains = 1,
	  parametersToSave = parameters.to.save,
	  nIter = 1000, nBurnin = 50000, nThin = 200,
	  engine = "WinBUGS", DIC=F)
	
	# compare the sample distributions for the alpha and tau
	# parameters but do not do so for the others; the 
	# posterior inferences are close but not exact,
	# due to the prior approximation here.
	# Adjust for the multiple test scripts being run.
	bugs.index <- c(1:4,60:61)
	bhlm.index <- (2:7)
	compareSampleDistributions(getSamples(bugs.samples, bugs.index), 
	  getSamples(bhlm.samples, bhlm.index), print.ks=F, 
	  pvalue.cutoff = (0.01 / 100))
}
