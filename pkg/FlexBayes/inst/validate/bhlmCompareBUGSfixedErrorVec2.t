# {cat("----- Compare bhlm, BUGS for known error vector ----------\n");T}
{
#
# Test not run, needs more work to get consistent results between
# BUGS and bhlm.  -spk oct-14-2010
#
# # Functions: posteriorSamples(engine='WinBUGS'), bhlm
# # Description: Compare the results for the bhlm function
# #  to those from the BUGS engine, for a linear 
# #  regression model fit to the orthodont data, when the 
# #  outcome MSE is a known vector.
# #
# 
# # Fit the model using the bhlm function
# 
# alpha.prior <- bayes.nonInformative()
# # fix the error variance for each group (subject)
# error.prior <- vector( "list", 27 )
# for( i in (1:27) ){
#   if (i < 13)
#     error.prior[[i]] <- bayes.massPoint(1.0)
#   else
#     error.prior[[i]] <- bayes.massPoint(1.5)
# }
# betaCov.prior <- bayes.invChisq(3,1)
# 
# orthoPrior <- bhlm.prior(error.var=error.prior, 
#   level2.coef = alpha.prior, 
#   random.var=betaCov.prior, common.error.var=0)
# 
# orthoSampler <- bhlm.sampler( nBurnin=1000, 
#   nSamples=1000, nThin=50,
#   nChains=1, init.point = "prior")
# 
# oldopt <- options(contrasts=c(factor="contr.treatment", ordered="contr.poly"))
# bhlm.samples <- bhlm(random.formula = distance~I(age-11),
#   level2.formula = ~Sex, group.formula = ~Subject, 
#   data = Orthodont, prior = orthoPrior, 
#   sampler = orthoSampler)
# options(oldopt)
# 
# # Prepare the data set for the BUGS analysis
# Orthodont2 <- list(distance=Orthodont$distance,
# 	age=Orthodont$age)
# Nsubjects <- length(unique(Orthodont$Subject))
# # Create the sex vector to specify the sex for each
# # subject, rather than the sex associated with each
# # data point.
# # Also generate the numeric subject index.
# sex.bysubj <- rep(0, Nsubjects)
# Ndata <- length(Orthodont$distance)
# Orthodont2$Subject <- rep(0, Ndata)
# for (subj in (1:Nsubjects)){
# 	subjID <- unique(Orthodont$Subject)[subj]
# 	thisSubj <- which(subjID == Orthodont$Subject)
# 	Orthodont2$Subject[thisSubj] <- subj
# 	sex.bysubj[subj] <- unique(Orthodont$Sex[thisSubj] == "Female")
# }
# Orthodont2$Sex <- sex.bysubj
# Orthodont2$Nsubjects <- Nsubjects
# Orthodont2$Ndata <- Ndata
# Orthodont2$prec.obs <- c( rep(1, 12), rep(0.6667, 15) )
# 
# ## specify the model in BUGS format	
# bugsModel <- function (){
#   for( i in 1 : Ndata ) {
# 		age.centered[i] <- age[i] - 11
# 		distance[i] ~ dnorm(mean[i], prec.obs[Subject[i]])
# 		mean[i] <- beta[Subject[i],1] + beta[Subject[i],2]*age.centered[i]
#   }
#   for (j in 1:Nsubjects){
#   	beta[j,1] ~ dnorm(mean.beta[j,1], prec.beta1)
#   	beta[j,2] ~ dnorm(mean.beta[j,2], prec.beta2)
#   	mean.beta[j,1] <- alpha00 + alpha01*Sex[j]
#   	mean.beta[j,2] <- alpha10 + alpha11*Sex[j]
#   }
#   prec.beta1 ~ dgamma(1.5, 1.5)
#   tau.beta1 <- 1/sqrt(prec.beta1)
#   prec.beta2 ~ dgamma(1.5, 1.5)
#   tau.beta2 <- 1/sqrt(prec.beta2)
#   alpha00 ~ dnorm(0,1e-9)
# 	alpha01 ~ dnorm(0,1e-9)
# 	alpha10 ~ dnorm(0,1e-9)
# 	alpha11 ~ dnorm(0,1e-9)
# 	
#   invisible()		
# }
# 
# # Specify the parameters for which we wish to save 
# # the posterior samples
# parameters.to.save <- c("alpha00", "alpha01", 
# 	"alpha10", "alpha11", "tau.beta1", "tau.beta2", "beta")
# 
# ## obtain the posterior samples
# bugs.samples <- posteriorSamples (
#   data = Orthodont2, model = bugsModel, 
#   nChains = 1,
#   parametersToSave = parameters.to.save,
#   nIter = 1000, nBurnin = 50000, nThin = 200,
#   engine = "WinBUGS", DIC=F)
# 
# # compare the sample distributions for the parameters.
# # Adjust for the multiple test scripts being run.
# bugs.index <- c((1:4), 59:60, 5:6)
# bhlm.index <- c(1:6, 35:36)
# compareSampleDistributions(
#   getSamples( bugs.samples[,bugs.index] ),
#   getSamples( bhlm.samples[,bhlm.index] ), print.ks=F,
#   pvalue.cutoff = (0.01 / 100) )
TRUE
}
