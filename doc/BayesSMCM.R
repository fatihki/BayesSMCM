## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(BayesSMCM)

## -----------------------------------------------------------------------------
# simulate data
n <- 300
b.true <- c( 0.4, 0.5, 0.1)
beta.true <- c(1, 0.2 )
baseline.hazard.rates <- 1
intervals <- c(0, 16)
prob.cov.X <- c(0.5)
prob.cov.Z <- c(0.5)
dat1 <- simSMCM(n, b.true, beta.true, baseline.hazard.rates, intervals, seed = 2026, cens.start = 30/365, cens.end = 32,  prob.cov.X, prob.cov.Z, same.cov = TRUE )
str(dat1)

## -----------------------------------------------------------------------------
# Bayesian settings 
nchains = 2; nIter   = 3000;   warmup = 1000; thin = 10

# prior for SMCMs
priorPar.smcm = list( r1 = 1, delta1 = 0.0001, r2 = 1, delta2 = 0.0001, a = 0.1, b = 0.1 )

# prior settings for SMCFMs
priorPar.smcfm = list( r1 = 1, delta1 = 0.0001, r2 = 1, delta2 = 0.0001, a = 0.1, b = 0.1, c = 0.1, d = 0.1 )

## -----------------------------------------------------------------------------
out.smcm.mcmc = fit.MCMC.SMCM ( data = dat1, hyperpar = priorPar.smcm, nchains, nIter, warmup, thin, mcmc.parallel = "parLapply", standardize = TRUE, probs = 0, save_loglik = 1, seed = 2025 )

## -----------------------------------------------------------------------------
print_smcm( out.smcm.mcmc, stan.model = FALSE, frailty = FALSE )

## -----------------------------------------------------------------------------
out.smcfm.mcmc = fit.MCMC.SMCFM ( data = dat1, hyperpar = priorPar.smcfm, nchains, nIter, warmup, thin, mcmc.parallel = "parLapply", standardize = TRUE, probs = 0, save_loglik = 1, seed = 2025 )

## -----------------------------------------------------------------------------
print_smcm( out.smcfm.mcmc, stan.model = FALSE, frailty = TRUE )

## -----------------------------------------------------------------------------
priorPar.smcm.rstan = list( sigma_beta = 1000, sigma_b = 1000, a = 0.1, b = 0.1)

## ----include=FALSE------------------------------------------------------------
out.smcm.rstan = fit.SMCM.RStan(data=dat1, hyperpar =  priorPar.smcm.rstan, nchains, nIter, warmup, thin, standardize = TRUE, probs = 0, save_loglik = 1, seed = 2025  )

## -----------------------------------------------------------------------------
print_smcm( out.smcm.rstan, stan.model = TRUE, frailty = FALSE)

## -----------------------------------------------------------------------------
priorPar.smcfm.rstan = list( sigma_beta = 1000, sigma_b = 1000,a = 0.1, b = 0.1,  c = 0.1, d = 0.1 )

## ----include=FALSE------------------------------------------------------------
out.smcfm.rstan =  fit.SMCFM.RStan(data=dat1, hyperpar =  priorPar.smcfm.rstan, nchains, nIter, warmup, thin, standardize = TRUE, probs = 0, save_loglik = 1, seed = 2025  )

## -----------------------------------------------------------------------------
print_smcm( out.smcfm.rstan, stan.model = TRUE, frailty = TRUE)

## ----include=FALSE------------------------------------------------------------
out.hsmcm.rstan = fit.HSMCM.RStan(data=dat1, hyperpar =  priorPar.smcm , nchains, nIter, warmup, thin, standardize = TRUE, probs = 0, save_loglik = 1, seed = 2025  )

## -----------------------------------------------------------------------------
print_smcm( out.hsmcm.rstan, stan.model = TRUE, frailty = FALSE)

## ----include=FALSE------------------------------------------------------------
out.hsmcfm.rstan = fit.HSMCFM.RStan(data=dat1, hyperpar = priorPar.smcfm, nchains, nIter, warmup, thin, standardize = TRUE, probs = 0, save_loglik = 1, seed = 2025  )

## -----------------------------------------------------------------------------
print_smcm( out.hsmcfm.rstan, stan.model = TRUE, frailty = TRUE)

