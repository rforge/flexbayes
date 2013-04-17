#########################################################################
##  
##  CGR (Cook, Gelman, and Rubin) Validation of FlexBayes 
##
##  
##

cgrFlexBayes <- function(modelName= "LinearReg", verbose = F, ...){
	
#########################################################################
##  
##  cgrFlexBayes
##
##  validates a specified model in FlexBayes using the method and code of 
##  Cook, Gelman, and Rubin (2006).
##  

	val <- switch (modelName, 
		LinearReg = validateLinearReg (verbose = verbose, ...),
		LogisticReg = validateLogisticReg (verbose = verbose, ...),
		LogisticTprior = validateLogisticTprior (verbose = verbose, ...),
		LogisticOverdisperse = validateLogisticOverdisperse (verbose = verbose, ...),
		PoissonGLM = validatePoissonGLM (verbose = verbose, ...),
		PoissonTprior = validatePoissonTprior (verbose = verbose, ...),
		PoissonOverdisperse = validatePoissonOverdisperse (verbose = verbose, ...),
		AEmodel = validateAEmodel (verbose = verbose, ...),
		# if modelName is not one of these
		stop("modelName must be one of: 'LinearReg', 'LogisticReg', 
     'LogisticTprior', 'LogisticOverdisperse',
		 'PoissonGLM', 'PoissonTprior', 'PoissonOverdisperse', 
		 or 'AEmodel'"))
	
	if (verbose)
		return(invisible(val))
	else  return(val)
}


rinvchisq <- function(n, df, sigmasq){

#########################################################################
##  
##  rinvchisq
##
##  generates from the scaled inverse chisquared distribution with degrees
##  of freedom df and squared scale sigmasq.  
##  

	alpha <- df/2
	beta <- alpha*sigmasq;
	draws <- 1/rgamma(n,alpha,beta)
	return(draws)
}

quant <- function (draws) 
{

#########################################################################
##  
##  quant
##
##  Calculates the empirical quantile of the first element of draws in 
##  the vector draws.  
##  

    n <- length(draws)
    rank.theta <- c(1:n)[order(draws) == 1] - 1
    quantile.theta <- (rank.theta + 0.5)/n
    return(quantile.theta)
}


bayes.validate <- function (generate.param, generate.param.inputs = NULL, generate.data, 
    generate.data.inputs = NULL, analyze.data, analyze.data.inputs = NULL, 
    n.rep = 20, n.batch = NULL, params.batch = NULL, print.reps = FALSE,
    quantilePlot.directory = NULL) 
{

#########################################################################
##  
##  bayes.validate
##
##  Runs the Cook, Gelman, and Rubin (2006) validation method.
##  

    if (is.null(generate.param.inputs)) 
        n.param <- length(generate.param())
    else n.param <- length(generate.param(generate.param.inputs))
    if (!is.null(n.batch)) {
        num.batches <- length(n.batch)
        if (sum(n.batch) != n.param) {
            print("Error:  Lengths of parameter batches don't sum to length of parameter vector")
            return()
        }
        if (!is.null(params.batch)) {
            if (length(params.batch) != length(n.batch)) {
                print("Error:  Must have same number of named parameters as batches")
                return()
            }
        }
        batch.ind <- rep(0, (num.batches + 1))
        for (i in 1:num.batches) batch.ind[(i + 1)] <- batch.ind[i] + 
            n.batch[i]
        plot.batch <- rep(1, n.batch[1])
        for (i in 2:num.batches) plot.batch <- c(plot.batch, 
            rep(i, n.batch[i]))
    }
    if (!is.null(n.batch)) 
        quantile.theta <- matrix(0, nrow = n.rep, ncol = (n.param + 
            num.batches))
    else quantile.theta <- matrix(0, nrow = n.rep, ncol = n.param)
    
    for (reps in 1:n.rep) {
        if (is.null(generate.param.inputs)) 
            theta.true <- generate.param()
        else theta.true <- generate.param(generate.param.inputs)
        if (is.null(generate.data.inputs)) 
            data.rep <- generate.data(theta.true)
        else data.rep <- generate.data(theta.true, generate.data.inputs)
        if (is.null(analyze.data.inputs)) 
            theta.draws <- analyze.data(data.rep, theta.true)
        else theta.draws <- analyze.data(data.rep, theta.true, 
            analyze.data.inputs)
        if (is.matrix(theta.draws) == FALSE) 
            theta.draws <- as.matrix(theta.draws)
        if (length(theta.true) != ncol(theta.draws)) {
            print("Error: Generate.param and analyze.data must output same number of parameters")
            print("True # params:")
						print(length(theta.true))
						print("# sampled params")
						print(ncol(theta.draws))
						
            return()
        }
        if (!is.null(n.batch)) {
            for (i in 1:num.batches) {
                if (n.batch[i] > 1) {
                  theta.draws <- cbind(theta.draws, apply(as.matrix(theta.draws[, 
                    (batch.ind[i] + 1):batch.ind[(i + 1)]]), 
                    1, mean))
                  theta.true <- c(theta.true, mean(theta.true[(batch.ind[i] + 
                    1):batch.ind[(i + 1)]]))
                }
                else {
                  theta.draws <- cbind(theta.draws, theta.draws[, 
                    (batch.ind[i] + 1)])
                  theta.true <- c(theta.true, theta.true[(batch.ind[i] + 
                    1)])
                }
            }
        }
        theta.draws <- rbind(theta.true, theta.draws)
        quantile.theta[reps, ] <- apply(theta.draws, 2, quant)
        if (print.reps == TRUE) 
            print(reps)
    }

		# quantile.theta has one column for each parameter plus one 
		# column for each batch, if n.batch is not null.
    quantile.trans1 <- (apply(quantile.theta, 2, qnorm))
    quantile.trans <- quantile.trans1^2
    q.trans <- apply(quantile.trans, 2, sum)
    pvals.chisq.right <- 1 - pchisq(q.trans, df = n.rep)
    pvals.chisq.left <- pchisq(q.trans, df = n.rep)
    if (is.null(n.batch)) {
        adj.min.pchisq.right <- n.param * min(pvals.chisq.right)
        adj.min.pchisq.left <- n.param * min(pvals.chisq.left)
    }
    else {  ## needs to be fixed
        z.batch <- z.stats[(n.param + 1):length(p.vals.chisq)]
    		# just take the p-values corresponding to the batches,
    		# and ignore those corresponding to the individual variables
        p.batch <- p.vals.chisq[(n.param + 1):length(p.vals.chisq)]
        adj.min.p.chisq <- num.batches * min(p.batch)
    }
    
    # now calculate the p-values from a normality test
   	pvals.normal <- KSnormalityTest (variables = quantile.trans1)
   	adj.min.pnormal <- n.param * min(pvals.normal)
    
    # plot histograms of the quantiles  
    if (!is.null(quantilePlot.directory))
    	plotQuantiles (quantiles = quantile.theta, 
  			directory = quantilePlot.directory)
  			
    # do a normality test on the transformed quantiles 
 
    if (is.null(n.batch)) {
        return(list(pvals.chisq.right = pvals.chisq.right, 
        pvals.chisq.left = pvals.chisq.left,
        pvals.normal = pvals.normal, 
        adj.min.pchisq.right = adj.min.pchisq.right, 
        adj.min.pchisq.left = adj.min.pchisq.left, 
        adj.min.pnormal = adj.min.pnormal))
    }
    else {  ## needs to be fixed
        if (length(z.batch) == n.param) 
            return(list(p.batch = p.batch, adj.min.p = adj.min.p))
        else return(list(p.vals = p.vals[1:n.param], p.batch = p.batch, 
            adj.min.p = adj.min.p))
    }
}


plotQuantiles <- function (quantiles, directory){

#########################################################################
##  
##  plotQuantiles
##
##  Plots a histogram of the quantiles for each parameter, which should 
##  appear uniform if posterior samples are being drawn from the correct
##  distribution.  Also transforms the quantiles via the normal inverse
##  CDF function \Phi^{-1}, so that the resulting values should be 
##  normally distributed, then creates a qq-plot of these values.
##  

  nParams <- dim(quantiles)[2]
  nPages <- ceiling(nParams/25)  # Only 25 plots fit comfortably
                                        # on a page
  quantile.trans <- apply(quantiles, 2, qnorm)

  for (page in (1:nPages)){
    filename <- paste (directory, "/uniformPlots", page, ".ps", sep = "")
    postscript(filename)
    par(mfrow = c(5, 5))
    firstParam <- (page-1)*25+1
    lastParam <- min(page*25,nParams)
    for (param in (firstParam:lastParam)){
      main = paste ("Parameter ", param, sep = "")
      hist(x = quantiles[,param], breaks =10, main = main,
         xlab = "quantile", ylab = "frequency")
    }
    dev.off()

    filename <- paste (directory, "/qqPlots", page, ".ps", sep = "")
    postscript(filename)
    par(mfrow = c(5, 5))
    firstParam <- (page-1)*25+1
    lastParam <- min(page*25,nParams)
    for (param in (firstParam:lastParam)){
      main = paste ("Parameter ", param, sep = "")
      qqnorm(x = quantile.trans[,param], main = main)
      qqline(x = quantile.trans[,param])
    }
    dev.off()
  }

}


KSnormalityTest <- function (variables){

#########################################################################
##  
##  KSnormalityTest
##
##  Tests whether the variables are distributed normally.
##  Each column is tested separately.
##  
	
	nParams <- dim(variables)[2]
	nDraws <- dim(variables)[1]
	
	ks.pvals <- rep(1,nParams)
	if (nDraws > 2){
		for (i in 1:nParams){
			if (is.R()){
				ks <- ks.test (x = variables[,i], y = "pnorm", mean = 0, sd = 1)
			} else {
				ks <- ks.gof(variables[,i], distribution = "normal", mean = 0, sd = 1)		
			}
			ks.pvals[i] <- ks$p.value
		}
	}
	return (ks.pvals)
}

printPvalues <- function(pvals, adj.min.pval, threshold){
 
  cat ("alpha-level: ", threshold, "\n")
  cat ("Unadjusted p-values for each of the parameters:\n")
  cat(signif(pvals, digits=3))
  cat ("\nBonferroni-adjusted minimum p-value:\n")
  cat(signif(adj.min.pval, digits=3))
  cat ("\n")
}

# pretty-print CGR test output
printCGRoutput <- function (test, threshold){
	failed.chisq.right <- (test$adj.min.pchisq.right < threshold)
	failed.chisq.left <- (test$adj.min.pchisq.left < threshold)
  failed.normal <- (test$adj.min.pnormal < threshold)

  cat("----------------- RESULTS -----------------------------------\n")	
	if (!failed.chisq.right && !failed.chisq.left && !failed.normal)
	  cat ("All tests PASSED\n")
  cat("----------------- Right-tail Chisquared Test ----------------\n")	
  if (failed.chisq.right){
  	cat ("Test FAILED\n")
  } else {
    cat ("Test PASSED\n")
  }
  printPvalues(test$pvals.chisq.right, test$adj.min.pchisq.right, 
		threshold)
  cat("----------------- Left-tail Chisquared Test -----------------\n")	
  if (failed.chisq.left){
  	cat ("Test FAILED\n")
  } else {
    cat ("Test PASSED\n")
  }
  printPvalues(test$pvals.chisq.left, test$adj.min.pchisq.left, 
		threshold)
  cat("----------------- Normality Test ----------------------------\n")	
  if (failed.normal){
  	cat ("Test FAILED\n")
  } else {
    cat ("Test PASSED\n")
  }
	printPvalues(test$pvals.normal, test$adj.min.pnormal, 
		threshold)
	
	return()
}


