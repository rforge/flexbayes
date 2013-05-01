summary.posterior <- function(object, ...)
{
  call <- object$call
  mcs <- summary(object$chains)
  DIC <- object$DIC
  ans <- list(call = call, mcs = mcs, DIC = DIC)
	oldClass(ans) <- "summary.blm"
	ans
}


print.summary.posterior <- function(x, ...)
{
  cat("\n  *** Posterior Distribution from the Bayesian Model ***\n")
  if(!is.null(x$call)){
    cat("\nCall:\n")
    print(x$call)
  }

  print(x$mcs, ...)

	if(!is.null(x$DIC)) {
		cat("3. DIC statistics:\n")
	  print(x$DIC, ...)
	}

  invisible(x)
}


print.posterior <-  function(x, ...){
	print(summary(x))
}


plot.posterior <- function(x, y, ...)
  plot(x$chains, ...)








