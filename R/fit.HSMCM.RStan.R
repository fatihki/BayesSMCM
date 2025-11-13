#' @title  Fit Bayesian Hierarchical Semi-parametric Mixture Cure Model (HSMCM) with RStan
#'
#' @description
#' This function fits the Bayesian HSMCM model using \pkg{RStan}, 
#' employing hierarchical priors analogous to the Bayesian lasso framework 
#' for both \eqn{\bm b} and \eqn{\bm \beta} regression coefficients as in 
#' Kızılaslan and Vitelli (2025).
#'
#' @name fit.HSMCM.RStan
#'  
#' @import rstan
#'
#' @param data a list containing observed data from \code{n} subjects,
#' \code{observed_time}, \code{delta}, \code{X}, \code{Z}. See details for more information.
#' @param hyperpar a list containing prior parameter values for the HSMCM model; among c("r1", "delta1", "r2", "delta2", "a", "b").
#' See details for more information.
#' @param nchains Number of MCMC chains.
#' @param nIter Number of MCMC iterations.
#' @param warmup Number of warmup iterations.
#' @param thin Thinning interval for MCMC samples.
#' @param standardize Logical indicating whether to standardize covariates.
#' @param probs A numeric vector of quantiles used to define the cut points of the piecewise exponential baseline hazard function. 
#' For example, use \code{probs = 0} when \eqn{J = 1}, \code{probs = 0.5} or any other quantile value when \eqn{J = 2}, 
#' and a vector such as \code{probs = c(0.35, 0.70)} when \eqn{J = 3}.
#' @param save_loglik Integer indicating whether to save log-likelihood values (0 = FALSE, 1 = TRUE).
#' @param seed Random seed for reproducibility.
#' 
#' @details
#' \tabular{ll}{
#' \code{observed_time} \tab a vector of \eqn{n} times to the event.\cr
#' \code{delta} \tab a vector of \eqn{n} censoring indicators for the event time (1=event occurred, 0=censored).\cr
#' \code{X} \tab covariate matrix for the incidence part of the model \eqn{\bm X}, \eqn{n} observations by \eqn{p_2} variables,
#' where \eqn{p_2} is the length of \eqn{\bm \beta}.\cr
#' \code{Z} \tab covariate matrix for the latency part of the model \eqn{\bm Z}, \eqn{n} observations by \eqn{p_1+1} variables,
#' where \eqn{p_1+1} is the length of \eqn{\bm \b}.\cr
#' \code{r1} \tab the shape parameter of the gamma prior for \eqn{\eta^2}.\cr
#' \code{delta1} \tab the rate parameter of the gamma prior for \eqn{\eta^2}.\cr
#' \code{r2} \tab the shape parameter of the gamma prior for \eqn{\eta^{*2}}.\cr 
#' \code{delta2} \tab the rate parameter of the gamma prior for \eqn{\eta^{*2}}.\cr
#' \code{a} \tab the shape parameter of the gamma prior for \eqn{\bm \lambda}.\cr 
#' \code{b} \tab the rate parameter of the gamma prior for \eqn{\bm \lambda}.\cr
#' }
#' 
#' @return A list containing:
#' \item{fit}{A stanfit object.}
#' \item{fit.results}{Summary of posterior estimation results, including posterior means, credible intervals, and 
#' convergence diagnostics for all model parameters.}
#' \item{b.chains}{Posterior samples for the regression coefficients \eqn{\bm b} across all MCMC chains.}
#' \item{beta.chains}{Posterior samples for the regression coefficients \eqn{\bm \beta} across all MCMC chains.}
#' \item{lambda.chains}{Posterior samples for the regression coefficients \eqn{\bm \lambda} across all MCMC chains.}
#' \item{Rhat}{Gelman–Rubin convergence statistic for assessing MCMC chain convergence.}
#' \item{loo.results}{Model fit statistics including DIC and LPML when \code{save_loglik = TRUE}.}
#' \item{loglik.matrix.chains}{Matrix of the log-likelihood of the model over MCMC samples when \code{save_loglik = TRUE}.}
#' \item{data}{Processed survival data used for model fitting.}
#' \item{priors}{Prior parameter settings used in the model.}
#' \item{mcmc.info}{nchains, nIter, warmup and thin are used in the model.}
#' 
#' @examples
#' # Fit HSMCM RStan model to data generated from Scenario 1 in our simulation study
#' n <- 300
#' b.true <- c( 0.4, 0.5, 0.1)
#' beta.true <- c(1, 0.2 )
#' baseline.hazard.rates <- 1
#' intervals <- c(0, 16)
#' prob.cov.X <- c(0.5)
#' prob.cov.Z <- c(0.5)
#' dat1 <- simSMCM(n, b.true, beta.true, baseline.hazard.rates, intervals, seed = 2025, 
#'   cens.start = 30/365, cens.end = 32,  prob.cov.X, prob.cov.Z, same.cov = TRUE )
#' str(dat1)
#' nchains = 2; nIter = 2500;   warmup = 500; thin = 10
#' priorPar = list( r1 = 1, delta1 = 0.0001, r2 = 1, delta2 = 0.0001, a = 0.1, b = 0.1 )
#' out.hsmcm.rstan =  fit.HSMCM.RStan(data = dat1, hyperpar =  priorPar, nchains, nIter, warmup, 
#'   thin, standardize = FALSE, probs = 0 , save_loglik = 1, seed = 2025  )
#' print_smcm( out.hsmcm.rstan, stan.model = TRUE, frailty = FALSE)
#' 
#' @references 
#' Kızılaslan, F., Vitelli, V. (2025). Bayesian Semiparametric Mixture Cure (Frailty) Model, arXiv \url{https://arxiv.org/???}.
#' 
#' @export
fit.HSMCM.RStan           <- function(data, hyperpar = list(), nchains = 2, nIter = 6000, warmup = 1000, thin = 20,
                                      standardize,  probs, save_loglik, seed = NULL  ){
  
  if (!is.null(seed)) set.seed(seed)
  seed.sampling          <- ifelse((!is.null(seed)), seed, NA)
  
  n                      <- length(data$observed_time)
  p1 <- ncol(data$X); p2 <- ncol(data$Z)
  
  survObj                <- list()
  survObj$t              <- data$observed_time
  survObj$di             <- data$delta
  if( standardize == TRUE){  
    survObj$X            <- data$X
    survObj$Z            <- data$Z
  }else{
    survObj$X            <- scale.dummy.matrix(data$X)
    survObj$Z            <- scale.dummy.matrix(data$Z) 
  }
  
  event.time.points		  <- sort( survObj$t[survObj$di == 1] )
  if( probs[1] == 0 ){
    s	                  <- c( 0,  2*max(survObj$t) - max( survObj$t[-which(survObj$t==max(survObj$t))] ) )
  }else{
    s	                  <- c( 0, as.numeric( quantile(event.time.points, probs = probs ) ),  2*max(survObj$t) - max( survObj$t[-which(survObj$t==max(survObj$t))] ) )
  }
  JJ                    <- length(s)-1
  
  priorPara             <- list()
  priorPara$s           <- s 
  priorPara$J           <- JJ
  priorPara$r1          <- hyperpar$r1 
  priorPara$delta1      <- hyperpar$delta1 
  priorPara$r2          <- hyperpar$r2
  priorPara$delta2      <- hyperpar$delta2
  priorPara$a           <- hyperpar$a 
  priorPara$b           <- hyperpar$b 
  
  data.stan.model       <- list(n=n, p1=p1, p2=p2, JJ=JJ, s=s, t=survObj$t, delta=survObj$di, X=survObj$X, Z=survObj$Z,
                                lambda_a0 = priorPara$a, lambda_b0 = priorPara$b, 
                                r1 = priorPara$r1, delta1 = priorPara$delta1, 
                                r2 = priorPara$r2, delta2 = priorPara$delta2,
                                save_loglik = save_loglik )
  # HSMCM.model           <- stan_model(file = system.file("stan", "HSMCM.stan", package = "BayesSMCM") )
  fit                   <- sampling( HSMCM.model, data = data.stan.model,
                                     iter = nIter, warmup = warmup, thin = thin, chains = nchains, cores = nchains, refresh =  nIter/2,
                                     seed = seed.sampling,
                                     control = list(adapt_delta=0.80, stepsize=0.1, max_treedepth=15) )
  out                   <- summary(fit)$summary
  Rhat                  <- out[, "Rhat"] 
  if ( check.stan.rhat(fit=fit,  rhat_threshold = 1.1) == FALSE ){
    fit                   <- sampling( HSMCM.model, data = data.stan.model,
                                       iter = nIter, warmup = warmup, thin = thin, chains = nchains, cores = nchains, refresh =  nIter/2,
                                       seed = seed.sampling,
                                       control = list(adapt_delta=0.75, stepsize=0.1, max_treedepth=15) )
    out                 <- summary(fit)$summary
    Rhat                <- out[, "Rhat"] 
  }
  b.mcmc.chains         <- extract(fit)$b
  beta.mcmc.chains      <- extract(fit)$beta
  lambda.mcmc.chains    <- extract(fit)$lambda
  
  if(save_loglik == 1){
    loglik.matrix.chains<- extract(fit)$log_lik
    loo.results         <- mcmc.dic.lpml( survObj, priorPara, loglik.matrix.chains, b.hat.mcmc = colMeans(b.mcmc.chains),
                                          beta.hat.mcmc = colMeans(beta.mcmc.chains) , lambda.hat.mcmc = colMeans(lambda.mcmc.chains), 
                                          stan.model = TRUE, frailty = FALSE )
  }
  fit.results           <- c( b.hat = colMeans(b.mcmc.chains),
                              beta.hat = colMeans(beta.mcmc.chains), 
                              lambda.hat =  colMeans(lambda.mcmc.chains),
                              pz.hat = mean( logit( as.matrix(survObj$Z), colMeans(b.mcmc.chains) ) ),
                              time = sum(get_elapsed_time(fit))/60 )
  
  ret <- list(  fit = fit,
                fit.results = fit.results,
                b.chains = b.mcmc.chains,
                beta.chains = beta.mcmc.chains,
                lambda.chains = lambda.mcmc.chains,
                Rhat = Rhat,
                loo.results = if( save_loglik ==1 ) loo.results else NULL,
                loglik.matrix.chains = if( save_loglik ==1 ) loglik.matrix.chains else NULL,
                data = survObj,
                priors = priorPara,
                mcmc.info = list( nchains = nchains, nIter = nIter, warmup = warmup, thin = thin )
  )
  class(ret) <- "HSMCM.RStan"
  
  return(ret)  
  
}