
# Summarize the distributions of the parameters 
# in a posterior object

summary.posterior <- function(object, maxVars = 30, digits = 4, ...)
{
	nParams <- nvar(object)
	nParams <- min(nParams, maxVars)
	
	paramMean <- rep(0, nParams)
	paramSd <- rep(0, nParams)
	probs <- c(0.025, 0.25, 0.50, 0.75, 0.975)
	paramQuantiles <- matrix(0, nrow = nParams, ncol = length(probs))

  for (i in (1:nParams)){
		#samples <- getSamples(x=x, param=i)
    samples <- as.vector(sapply(object, function(u, idx) u[, idx], idx = i))
    paramMean[i] <- mean(samples)
		paramSd[i] <- sqrt(var(samples))
		paramQuantiles[i,] <- quantile(samples, probs = probs)
	}
	
	paramSummary <- data.frame(paramMean, paramSd, paramQuantiles)
	paramSummary <- signif(paramSummary, digits)
	if(!is.null(varnames(object)))
	  rownames(paramSummary) <- varnames(object)[1:nParams]
	probNames <- paste(as.character(100*probs), "%")
	colnames(paramSummary) <- c("Mean", "S.D.", probNames)
	ans <- list(paramSummary = paramSummary)

  if(!is.null(object$DIC))
		ans$DIC <- object$DIC

  ans$nchain <- nchain(object)
	ans$start <- start(object)
	ans$end <- end(object)
	ans$thin <- thin(object)
	ans$niter <- niter(object)
  ans$call <- object$call

  oldClass(ans) <- "summary.posterior"
	ans
}

print.summary.posterior <- function(x, ...){
  cat("*** Posterior Distribution from the Bayesian Model ***\n")
  if (!is.null(x$call)){
    cat("Call:  \n")
    print(x$call)
    cat("\n\n")
  }

	cat ("# of Chains: ", x$nchain, "\n")
	cat ("Starting Iteration: ", x$start, "\n")
	cat ("Ending Iteration: ", x$end, "\n")
	cat ("Thinning: ", x$thin, "\n")
	cat ("# of Samples: ", x$niter, "\n\n")
	
	cat("1. Summary statistics:\n\n")
	print(x$paramSummary[,1:2])
	cat("\n2. Quantiles:\n\n")
	print(x$paramSummary[,3:7])
	if (!is.null(x$DIC)){
		cat("\n3. DIC statistics:\n\n")
	  print(x$DIC)
	}
	return()
}

# The print function for objects of class posterior
print.posterior <- function(x, ...)
	print(summary(x))


