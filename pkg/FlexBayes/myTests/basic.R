#"norm:nonInformative:nonInformative"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

(x <- blm(MSFT - t90 ~ market - t90, data = largecap.ts))


#"normal:nonInformative:invChisq"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

returns <- microcap.ts[49:60, "GAIT"]
pr.var <- fbprior("invChisq", df = 7.83, sigma0.sq = 0.013)
my.prior <- blm.prior(priorSigma = pr.var)

(x <- blm(GAIT ~ 1, data = returns, prior = my.prior))


#"normal:normal:nonInformative"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

pr.coef <- fbprior("norm", mean = c(0, 1), S = diag(c(0.01^2, 0.5^2)))
my.prior <- blm.prior(priorBeta = pr.coef)
(x <- blm(MSFT - t90 ~ market - t90, data = largecap.ts, prior = my.prior))


#"normal:normal:invChisq"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

pr.coef <- fbprior("norm", mean = c(0, 1), S = diag(c(0.01^2, 0.5^2)))
pr.var <- fbprior("invChisq", df = 5, sigma0.sq = 0.007)
my.prior <- blm.prior(priorBeta = pr.coef, priorSigma = pr.var)
(x <- blm(MSFT - t90 ~ market - t90, data = largecap.ts, prior = my.prior))


#"normal:t:nonInformative"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

pr.coef <- fbprior("t", mean = c(0, 1), S = diag(c(0.01^2, 0.5^2)), df = 4)
my.prior <- blm.prior(priorBeta = pr.coef)
(x <- blm(MSFT - t90 ~ market - t90, data = largecap.ts, prior = my.prior))


#"normal:t:invChisq"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

pr.coef <- fbprior("t", mean = c(0, 1), S = diag(c(0.01^2, 0.5^2)), df = 4)
pr.var <- fbprior("invChisq", df = 5, sigma0.sq = 0.007)
my.prior <- blm.prior(priorBeta = pr.coef, priorSigma = pr.var)
(x <- blm(MSFT - t90 ~ market - t90, data = largecap.ts, prior = my.prior))


#"normal:normalMixture:nonInformative"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

returns <- smallcap.ts[49:60, "KRON"]

pr.coef <- fbprior("normmix",
                   mean = list(-0.1, 0.30),
                   S = list(0.01, 0.01),
                   w = c(0.5, 0.5))

my.prior <- blm.prior(priorBeta = pr.coef)
(x <- blm(KRON ~ 1, data = returns, prior = my.prior))

mu <- rep(1, nrow(returns))
(x <- blm(KRON ~ mu - 1, data = returns, prior = my.prior))


#"normal:normalMixture:invChisq"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

returns <- smallcap.ts[49:60, "KRON"]

pr.coef <- fbprior("normmix",
                   mean = list(0.02, 0.08),
                   S = list(0.02, 0.02),
                   w = c(0.65, 0.35))

pr.var <- fbprior("invChisq", df = 7.83, sigma0.sq = 0.013)

my.prior <- blm.prior(priorBeta = pr.coef, priorSigma = pr.var)
(x <- blm(KRON ~ 1, data = returns, prior = my.prior))


#"normal:tMixture:nonInformative"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

returns <- smallcap.ts[49:60, "KRON"]

pr.coef <- fbprior("tmix",
                   mean = list(0.02, 0.08),
                   S = list(0.02, 0.02),
                   w = c(0.65, 0.35),
                   df = c(5, 5))


my.prior <- blm.prior(priorBeta = pr.coef)
(x <- blm(KRON ~ 1, data = returns, prior = my.prior))


#"normal:tMixture:invChisq"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

pr.coef <- fbprior("tmix",
                   mean = list(0.02, 0.08),
                   S = list(0.02, 0.02),
                   w = c(0.65, 0.35),
                   df = c(5, 5))

pr.var <- fbprior("invChisq", df = 7.83, sigma0.sq = 0.013)

my.prior <- blm.prior(priorBeta = pr.coef, priorSigma = pr.var)

(x <- blm(MSFT - t90 ~ market - t90, prior = my.prior, data = largecap.ts))


#"t:nonInformative:nonInformative"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

my.lik <- blm.likelihood("t", df = 4)

(x <- blm(MSFT - t90 ~ market - t90, likelihood = my.lik, data = largecap.ts))


#"t:nonInformative:invChisq"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

returns <- microcap.ts[49:60, "GAIT"]
pr.var <- fbprior("invChisq", df = 7.83, sigma0.sq = 0.013)
my.prior <- blm.prior(priorSigma = pr.var)

my.lik <- blm.likelihood("t", df = 4)

(x <- blm(GAIT ~ 1, data = returns, prior = my.prior, likelihood = my.lik))


#"t:normal:nonInformative"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

pr.coef <- fbprior("norm", mean = c(0, 1), S = diag(c(0.01^2, 0.5^2)))
my.prior <- blm.prior(priorBeta = pr.coef)

my.lik <- blm.likelihood("t", df = 4)

(x <- blm(MSFT - t90 ~ market - t90, data = largecap.ts, prior = my.prior, likelihood = my.lik))


#"t:normal:invChisq"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

pr.coef <- fbprior("norm", mean = c(0, 1), S = diag(c(0.01^2, 0.5^2)))
pr.var <- fbprior("invChisq", df = 5, sigma0.sq = 0.007)
my.prior <- blm.prior(priorBeta = pr.coef, priorSigma = pr.var)

my.lik <- blm.likelihood("t", df = 4)

(x <- blm(MSFT - t90 ~ market - t90, data = largecap.ts, prior = my.prior, likelihood = my.lik))


#"t:t:nonInformative"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

pr.coef <- fbprior("t", mean = c(0, 1), S = diag(c(0.01^2, 0.5^2)), df = 4)
my.prior <- blm.prior(priorBeta = pr.coef)

my.lik <- blm.likelihood("t", df = 4)

(x <- blm(MSFT - t90 ~ market - t90, data = largecap.ts, prior = my.prior, likelihood = my.lik))


#"t:t:invChisq"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

pr.coef <- fbprior("t", mean = c(0, 1), S = diag(c(0.01^2, 0.5^2)), df = 4)
pr.var <- fbprior("invChisq", df = 5, sigma0.sq = 0.007)
my.prior <- blm.prior(priorBeta = pr.coef, priorSigma = pr.var)

my.lik <- blm.likelihood("t", df = 4)

(x <- blm(MSFT - t90 ~ market - t90, data = largecap.ts, prior = my.prior, likelihood = my.lik))


#"t:normalMixture:nonInformative"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

returns <- smallcap.ts[49:60, "KRON"]

pr.coef <- fbprior("normmix",
mean = list(0.02, 0.08),
S = list(0.02, 0.02),
w = c(0.65, 0.35))

my.prior <- blm.prior(priorBeta = pr.coef)

my.lik <- blm.likelihood("t", df = 4)

(x <- blm(KRON ~ 1, data = returns, prior = my.prior, likelihood = my.lik))


#"t:normalMixture:invChisq"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

returns <- smallcap.ts[49:60, "KRON"]

pr.coef <- fbprior("normmix",
mean = list(0.02, 0.08),
S = list(0.02, 0.02),
w = c(0.65, 0.35))

pr.var <- fbprior("invChisq", df = 7.83, sigma0.sq = 0.013)

my.prior <- blm.prior(priorBeta = pr.coef, priorSigma = pr.var)

my.lik <- blm.likelihood("t", df = 4)

(x <- blm(KRON ~ 1, data = returns, prior = my.prior, likelihood = my.lik))


#"t:tMixture:nonInformative"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

returns <- smallcap.ts[49:60, "KRON"]

pr.coef <- fbprior("tmix",
mean = list(0.02, 0.08),
S = list(0.02, 0.02),
w = c(0.65, 0.35),
df = c(5, 5))

my.prior <- blm.prior(priorBeta = pr.coef)

my.lik <- blm.likelihood("t", df = 4)

(x <- blm(KRON ~ 1, data = returns, prior = my.prior, likelihood = my.lik))


#"t:tMixture:invChisq"

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)
library(mpo)

pr.coef <- fbprior("tmix",
mean = list(0.02, 0.08),
S = list(0.02, 0.02),
w = c(0.65, 0.35),
df = c(5, 5))

pr.var <- fbprior("invChisq", df = 7.83, sigma0.sq = 0.013)

my.prior <- blm.prior(priorBeta = pr.coef, priorSigma = pr.var)

my.lik <- blm.likelihood("t", df = 4)

(x <- blm(MSFT - t90 ~ market - t90, prior = my.prior, data = largecap.ts, likelihood = my.lik))


################################################################################

rm(list = ls())
#library(FlexBayes, lib.loc = "testlib")
library(FlexBayes)

pr.coef <- fbprior("norm", mean = c(-40, 0, 1, 0), S = diag(c(4, 0.25, 0.25, 0.25)))
my.prior <- blm.prior(priorBeta = pr.coef)
(object <- blm(stack.loss ~ ., data = stackloss, prior = my.prior))




