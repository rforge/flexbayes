coef.blm <- function(object, LOC = mean, agg = TRUE, ...)
{
  chains <- object$coefficients
  nChains <- nchain(chains)
  coef.names <- dimnames(chains[[1]])[[2]]
  coefs <- array(NA, c(nChains, length(coef.names)))
  dimnames(coefs) <- list(paste("Chain ", 1:nChains, ":", sep = ""), coef.names)
  for(chain in 1:nChains)
    coefs[chain, ] <- apply(chains[[chain]], 2, LOC)

  if(agg) colMeans(coefs) else coefs
}


print.blm <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  if(length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else
    cat("No coefficients\n")

  cat("\n")
  invisible(x)
}


fitted.blm <- function(object, LOC = mean, ...)
{
  cf <- coef(object, LOC = LOC)
  X <- model.matrix(object)
  drop(X %*% cf)
}


model.frame.blm <- function(formula, ...)
  formula$model


model.matrix.blm <- function(object, ...)
  model.matrix.default(object = object$terms, data = object$model,
                       contrasts.arg = object$contrasts)


residuals.blm <- function(object, LOC = mean, ...)
{
  Y <- model.extract(model.frame(object), "response")
  drop(Y - fitted(object, LOC = LOC))
}


summary.blm <- function(object, LOC = mean, DISP = sd, ...)
{
    cf <- coef(object, LOC = LOC)
    coefficients <- matrix(NA, length(cf), 4)
    coefficients[, 1] <- cf
    coefficients[, 2] <- coef(object, LOC = DISP)
    dimnames(coefficients) <- list(names(cf), c("Estimate", "Std. Error", "", ""))
    
    sigma <- mean(sapply(object$sigma, DISP))
    
    res <- residuals(object, LOC = LOC)
    
    ans <- list(call = object$call, terms = object$terms, residuals = res,
    coefficients = coefficients, sigma = sigma)
    oldClass(ans) <- "summary.blm"
    ans
}


print.summary.blm <- function(x, digits = max(3L, getOption("digits") - 3L),
                              signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  resid <- x$residuals
  cat("Residuals:\n")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- structure(zapsmall(quantile(resid), digits + 1L), names = nam)
  print(rq, digits = digits, ...)

  cat("\nCoefficients:\n")
  coefs <- x$coefficients
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)

  cat("\nResidual standard error:", format(signif(x$sigma, digits)))
  cat("\n\n")

  invisible(x)
}














