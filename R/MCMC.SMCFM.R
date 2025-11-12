#' @title Function to MCMC algorithm for Semi-parametric Mixture Cure Frailty Model (SMCFM)
#'
#' @description
#' This function obtains posterior samples for all model parameters using an MCMC algorithm 
#' for a single chain when applying the \code{MCMC.SMCFM} framework  in Kızılaslan and Vitelli (2025) to the data. 
#'
#' @name MCMC.SMCFM
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar flush.console
#'
#' @param survObj a list containing observed data from \code{n} subjects;
#' \code{t}, \code{di}, \code{X}, \code{Z}. See details for more information.
#' @param priorPara a list containing prior parameter values for the SMCM model; among c("r1", "delta1", "r2", "delta2", "a", "b", "c", "d", "s", "J").
#' See details for more information.
#' @param initial Initial value of parameters.
#' @param nIter Number of MCMC iterations.
#' @param save_loglik Integer indicating whether to save log-likelihood values (0 = FALSE, 1 = TRUE).
#'
#' @details
#' The input data and prior parameters include the following:
#'  \tabular{ll}{
#'   \code{t}\tab A vector of \eqn{n} times to the event.\cr
#'   \code{di}\tab A vector of \eqn{n} censoring indicators (1 = event occurred, 0 = censored).\cr
#'   \code{X} \tab covariate matrix for the incidence part of the model \eqn{\bm X}, \eqn{n} observations by \eqn{p_2} variables,
#' where \eqn{p_2} is the length of \eqn{\bm \beta}.\cr
#'   \code{Z} \tab covariate matrix for the latency part of the model \eqn{\bm Z}, \eqn{n} observations by \eqn{p_1+1} variables,
#' where \eqn{p_1+1} is the length of \eqn{\bm \b}.\cr
#'   \code{r1}\tab Shape parameter of the gamma prior for \eqn{\eta^2}.\cr
#'   \code{delta1}\tab Rate parameter of the gamma prior for \eqn{\eta^2}.\cr
#'   \code{r2}\tab Shape parameter of the gamma prior for \eqn{\eta^{*2}}.\cr
#'   \code{delta2}\tab Rate parameter of the gamma prior for \eqn{\eta^{*2}}.\cr
#'   \code{a}\tab Shape parameter of the gamma prior for \eqn{\bm \lambda}.\cr
#'   \code{b}\tab Rate parameter of the gamma prior for \eqn{\bm \lambda}.\cr
#'   \code{c}\tab Shape parameter of the gamma prior for the frailty parameter \eqn{\theta}.\cr
#'   \code{d}\tab Rate parameter of the gamma prior for the frailty parameter \eqn{\theta}.\cr
#'   \code{s}\tab Interval boundaries for the piecewise exponential baseline hazard.\cr
#'   \code{J}\tab Total number of intervals for the piecewise exponential baseline hazard.\cr
#'   }
#' 
#' @return A list containing:
#'    \item{beta.p}{MCMC samples for the regression coefficients \eqn{\bm \beta}.}
#'    \item{b.p}{MCMC samples for the regression coefficients \eqn{\bm b}.}
#'    \item{tauSq.p}{posterior samples for \eqn{\tau^2}.}
#'    \item{tauStarSq.p}{posterior samples for \eqn{\tau^{*2}}.}
#'    \item{lambda.p}{posterior samples for \eqn{\bm \lambda}.}
#'    \item{theta.p}{posterior samples for the frailty parameter \eqn{\theta}.}
#'    \item{mcmcOutcome}{The list containing posterior samples for the remaining model parameters.}
#'    \item{t}{Processed times to the event data.}
#'    \item{di}{Processed censoring indicators for the event time.}
#'    \item{loglik.matrix}{Matrix of the log-likelihood of the model over MCMC samples when \code{save_loglik = TRUE}.}
#' 
#' @references 
#' Kızılaslan, F., Vitelli, V. (2025). Bayesian Semiparametric Mixture Cure (Frailty) Model, arXiv \url{https://arxiv.org/???}.
#' 
#' @export
MCMC.SMCFM                  <- function(survObj, priorPara, initial, nIter = 6000, save_loglik = 0 ){
  
  r1		       	            <- priorPara$r1
  delta1			              <- priorPara$delta1
  r2		       	            <- priorPara$r2
  delta2			              <- priorPara$delta2
  a                         <- priorPara$a
  b                         <- priorPara$b
  c                         <- priorPara$c
  d                         <- priorPara$d
  s							            <- priorPara$s
  J							            <- priorPara$J	
  
  ini	                      <- initial
  beta.ini			            <- ini$beta.ini
  b.ini			                <- ini$b.ini
  sigmaSq				            <- ini$sigmaSq
  sigmaStarSq	              <- ini$sigmaStarSq
  etaSq			                <- ini$etaSq
  etaStarSq			            <- ini$etaStarSq
  tauSq			  	            <- ini$tauSq
  tauStarSq			            <- ini$tauStarSq
  lambda.ini		            <- ini$lambda.ini	
  theta.ini                 <- ini$theta.in
  
  mcmcPara                  <- list()
  mcmcPara$beta.prop.mean 	<- beta.ini
  mcmcPara$b.prop.mean	    <- b.ini
  mcmcPara$lambda.prop.mean	<- lambda.ini
  
  ini$sd.beta	              <- sqrt(sigmaStarSq*tauStarSq) 
  ini$sd.b	                <- sqrt(sigmaSq*tauSq)
  ini$Xbeta                 <- as.vector(survObj$X %*% beta.ini)
  ini$Zb	                  <- as.vector(survObj$Z %*% b.ini)
  
  sd.beta	                  <- ini$sd.beta
  sd.b	                    <- ini$sd.b	 
  
  ## for posterior samples
  mcmcOutcome	              <- list()
  mcmcOutcome$initial			  <- initial
  mcmcOutcome$priorPara		  <- priorPara
  
  beta.p						        <- beta.ini
  b.p						            <- b.ini
  lambda.p					        <- lambda.ini
  theta.p                   <- theta.ini
  tauSq.p						        <- tauSq
  tauStarSq.p						    <- tauStarSq
  mcmcOutcome$sigmaSq.p	    <- sigmaSq
  mcmcOutcome$sigmaStarSq.p	<- sigmaStarSq
  mcmcOutcome$etaSq.p		    <- etaSq
  mcmcOutcome$etaStarSq.p	  <- etaStarSq
  mcmcOutcome$accept.beta	  <- c( rep(0, length(beta.ini)) )
  mcmcOutcome$accept.b		  <- c( rep(0, length(b.ini)) )
  mcmcOutcome$accept.lambda	<- c( rep(0, length(lambda.ini)) )
  mcmcOutcome$accept.theta	<- c( rep(0, length(theta.ini)) )
  if (save_loglik == 1) { 
    loglik.p                  <- logLik_SMCFM(survObj, priorPara, b.coef = b.p, beta.coef = beta.p, lambda.coef = lambda.p,  theta.coef = theta.p ) } # for initial values
  
  pb  <- txtProgressBar(min = 0, max = nIter, style = 3)
  
  # MCMC sampling
  for(M in 1:nIter){

    # Updating regression parameters b_k:
    sampleRP.b	          <- UpdateRP.b.SMCFM(survObj, priorPara, mcmcPara, ini)
    b.ini	                <- ini$b.ini	      <- sampleRP.b$b.ini
    Zb		                <- ini$Zb	          <- sampleRP.b$Zb
    mcmcOutcome$accept.b	<- mcmcOutcome$accept.b + sampleRP.b$accept
    
    # Updating regression parameters beta_k:
    sampleRP.beta	          <- UpdateRP.beta.SMCFM(survObj, priorPara, mcmcPara, ini)
    beta.ini	              <- ini$beta.ini	   <- sampleRP.beta$beta.ini
    Xbeta		                <- ini$Xbeta	     <- sampleRP.beta$Xbeta
    mcmcOutcome$accept.beta	<- mcmcOutcome$accept.beta + sampleRP.beta$accept
    
    # Updating baseline hazard paramaters lambda_k:
    sampleBH.lambda	         <- UpdateBH.lambda.SMCFM(survObj, priorPara, mcmcPara, ini)
    lambda.ini	             <- ini$lambda.ini	<- sampleBH.lambda$lambda.ini
    mcmcOutcome$accept.lambda	<- mcmcOutcome$accept.lambda + sampleBH.lambda$accept
    
    # Updating frailty parameter theta:
    sample.theta	           <- Update.theta.SMCFM(survObj, priorPara, mcmcPara, ini)
    theta.ini	               <- ini$theta.ini	  <- sample.theta$theta.ini
    mcmcOutcome$accept.theta <- mcmcOutcome$accept.theta + sample.theta$accept
    
    # Updating 1/tauSq
    tauSq	    <- ini$tauSq	        <- UpdateTau(survObj, priorPara, mcmcPara, ini)
    
    # Updating 1/tauStarSq
    tauStarSq	<- ini$tauStarSq	    <- UpdateTauStar(survObj, priorPara, mcmcPara, ini)
    
    
    # Updating sigmaSq
    sigmaSq	   <- ini$sigmaSq	      <- UpdateSigma(survObj, priorPara, mcmcPara, ini)
    
    # Updating sigmaStarSq
    sigmaStarSq	<- ini$sigmaStarSq	<- UpdateSigmaStar(survObj, priorPara, mcmcPara, ini)
    
    
    # Updating etaSq		
    etaSq	      <- ini$etaSq	      <- UpdateEta(survObj, priorPara, mcmcPara, ini)
    
    # Updating etaStarSq		
    etaStarSq	  <- ini$etaStarSq	  <- UpdateEtaStar(survObj, priorPara, mcmcPara, ini)
    
    ini$sd.beta <- sd.beta	        <- sqrt(sigmaStarSq*tauStarSq) 
    ini$sd.b	  <- sd.b             <- sqrt(sigmaSq*tauSq)
    
    if (save_loglik == 1) { 
      loglik.ini                <- logLik_SMCFM(survObj, priorPara, b.coef = b.ini, beta.coef = beta.ini, lambda.coef = lambda.ini, theta.coef = theta.ini ) }
    
    ###### storing posterior samples
    beta.p						        <- rbind(beta.p, beta.ini, deparse.level = 0)
    b.p						            <- rbind(b.p, b.ini, deparse.level = 0)
    lambda.p				          <- rbind(lambda.p, lambda.ini, deparse.level = 0)
    theta.p				            <- rbind(theta.p, theta.ini, deparse.level = 0)
    if (save_loglik == 1) {  
      loglik.p                <- rbind(loglik.p, loglik.ini, deparse.level = 0)  
    } 
    
    tauSq.p						        <- rbind(tauSq.p, tauSq, deparse.level = 0)
    tauStarSq.p						    <- rbind(tauStarSq.p, tauStarSq, deparse.level = 0)
    mcmcOutcome$sigmaSq.p	    <- c(mcmcOutcome$sigmaSq.p, sigmaSq)
    mcmcOutcome$sigmaStarSq.p	<- c(mcmcOutcome$sigmaStarSq.p, sigmaStarSq)
    mcmcOutcome$etaSq.p		    <- c(mcmcOutcome$etaSq.p, etaSq)
    mcmcOutcome$etaStarSq.p		<- c(mcmcOutcome$etaStarSq.p, etaStarSq)
    mcmcOutcome$ini				    <- ini
    
    # update progress bar every 10% of iterations
    if (M %% (nIter / 10) == 0 || M == nIter) {
      setTxtProgressBar(pb, M)
      flush.console()
    }
    
  } 
  
  close(pb)
  
  if (save_loglik == 1) { 
    ret      <- list(beta.p = beta.p, b.p = b.p, tauSq.p = tauSq.p,  tauStarSq.p = tauStarSq.p, lambda.p = lambda.p, theta.p = theta.p,  mcmcOutcome = mcmcOutcome, t = survObj$t, di = survObj$di, loglik.matrix = loglik.p )
  }else{
    ret      <- list(beta.p = beta.p, b.p = b.p, tauSq.p = tauSq.p,  tauStarSq.p = tauStarSq.p, lambda.p = lambda.p, theta.p = theta.p,  mcmcOutcome = mcmcOutcome, t = survObj$t, di = survObj$di ) 
  }
  
  return(ret)
} 
#'
UpdateRP.b.SMCFM <- function(survObj, priorPara, mcmcPara, ini){
  
  Z              <- survObj$Z
  X              <- survObj$X
  y              <- survObj$t
  delta          <- survObj$di
  J              <- priorPara$J
  v              <- priorPara$v
  s              <- priorPara$s
  b.ini          <- ini$b.ini	 
  beta.ini       <- ini$beta.ini
  lambda.ini     <- ini$lambda.ini
  theta.ini      <- ini$theta.ini
  sigmaSq        <- ini$sigmaSq	
  tauSq          <- ini$tauSq	
  sd.b           <- ini$sd.b	
  n              <- nrow(Z)
  p              <- length(b.ini)
  accept         <- rep(0, p)
  Zb             <- Z %*% b.ini
  Xbeta          <- X %*% beta.ini
  
  H0_h0                   <- cumulative_hazard(y, interval_bounds=s, lambda=lambda.ini)
  H0                      <- H0_h0$H0 
  h0                      <- H0_h0$h0 
  exp.Xbeta               <- exp(Xbeta)
  H1                      <- H0 * exp.Xbeta 
  
  for (k in 1:p) {
    
    exp.Zb                  <- exp(Zb)  
    
    term1                   <- Zb * delta 
    term2                   <- log1p(exp.Zb)
    H2                      <- 1 + (H1 / theta.ini)  
    term3                   <- log1p(exp.Zb * H2^(-theta.ini) )* (1-delta)  
    loglik.ini              <- sum( term1 - term2 + term3 )                     # log-likelihood for b_k at initial value
    
    # proposal value
    b.prop.mean             <- b.ini[k] 
    b.prop.variance         <- 1^2 
    b.prop                  <- b.ini
    b.prop[k]	              <- mean( rnorm(5, mean = b.prop.mean, sd = sqrt(b.prop.variance) ) )
    # Calculating acceptance probability
    Zb.prop	                <- Zb - Z[,k] * b.ini[k] + Z[,k] * b.prop[k]        # updating Zb with new b[k]
    exp.Zb.prop	            <- exp(Zb.prop)
    
    term1.prop              <- Zb.prop * delta
    term2.prop              <- log1p(exp.Zb.prop)
    term3.prop              <- log1p(exp.Zb.prop * (H2^(-theta.ini) ) )* (1-delta)
    loglik.prop             <- sum( term1.prop - term2.prop + term3.prop  )     # log-likelihood at proposal value of b[k]
    
    b.prop.mean.ini         <- b.prop[k] 
    b.prop.variance.ini     <- 1^2 
    #log-prior density at b.ini[k] and b.prop[k]
    logprior.ini            <- dnorm(b.ini[k],  mean = 0 , sd = sd.b[k], log = TRUE) 
    logprior.prop           <- dnorm(b.prop[k], mean = 0 , sd = sd.b[k], log = TRUE)
    # log-proposal density at b.ini[k] and b.prop[k]
    logprop.ini             <- dnorm(b.ini[k], mean = b.prop.mean, sd = sqrt(b.prop.variance), log = TRUE)
    logprop.prop            <- dnorm(b.prop[k], mean = b.prop.mean.ini, sd = sqrt(b.prop.variance.ini), log = TRUE)
    logR                    <- loglik.prop - loglik.ini + logprior.prop - logprior.ini + logprop.ini - logprop.prop
    
    u = log(runif(1)) < logR
    
    if(u == 1){
      b.ini[k]  <- b.prop[k]
      Zb	      <- Zb.prop
    }
    accept[k]   <- accept[k] + u
  } 
  list(b.ini = b.ini, accept = accept, Zb = Zb)
}
#'
UpdateRP.beta.SMCFM <- function(survObj, priorPara, mcmcPara, ini){
  
  X                <- survObj$X
  Z                <- survObj$Z
  y                <- survObj$t
  delta            <- survObj$di
  J                <- priorPara$J
  v                <- priorPara$v
  s                <- priorPara$s
  b.ini            <- ini$b.ini
  beta.ini         <- ini$beta.ini	 
  sigmaStarSq      <- ini$sigmaStarSq	
  tauStarSq        <- ini$tauStarSq
  lambda.ini       <- ini$lambda.ini
  theta.ini        <- ini$theta.ini
  sd.beta          <- ini$sd.beta	
  
  n                <- nrow(X)
  p                <- length(beta.ini)
  accept           <- rep(0, p)
  Xbeta            <- X %*% beta.ini 
  Zb               <- Z %*% b.ini 
  exp.Zb           <- exp(Zb)
  H0_h0            <- cumulative_hazard(y, interval_bounds=s, lambda=lambda.ini)
  H0               <- H0_h0$H0 
  h0               <- H0_h0$h0 
  
  for (k in 1:p) {
    
    exp.Xbeta                    <- exp(Xbeta)
    H1                           <- H0 * exp.Xbeta
    H20                          <- H1 / theta.ini  
    H2                           <- 1 + H20         
    
    term1                        <- (Xbeta - (theta.ini+1) * log1p(H20) ) * delta 
    term2                        <- log1p(exp.Zb * H2^(-theta.ini)) * (1-delta)  
    loglik.ini                   <- sum(term1 + term2) 
    
    beta.prop.mean               <- beta.ini[k]
    beta.prop.variance           <- 1^2
    beta.prop                    <- beta.ini
    beta.prop[k]	               <- mean( rnorm(5, mean = beta.prop.mean, sd = sqrt(beta.prop.variance) ) )
    # Calculating acceptance probability
    Xbeta.prop	                 <- Xbeta - X[,k] * beta.ini[k] + X[,k] * beta.prop[k] 
    exp.Xbeta.prop	             <- exp(Xbeta.prop)
    
    H1.prop                      <- H0 * exp.Xbeta.prop
    H20.prop                     <- H1.prop / theta.ini
    H2.prop                      <- 1 + H20.prop  
    term1.prop                   <- (Xbeta.prop - (theta.ini+1) * log1p(H20.prop) ) * delta
    term2.prop                   <- log1p(exp.Zb * H2.prop^(-theta.ini)) * (1-delta)
    loglik.prop                  <- sum(term1.prop + term2.prop) 
    beta.prop.mean.ini           <- beta.prop[k]
    beta.prop.variance.ini       <- 1^2
    #log-prior density at beta.ini[k] and beta.prop[k]
    logprior.ini                 <- dnorm(beta.ini[k],  mean = 0, sd = sd.beta[k], log = TRUE) 
    logprior.prop                <- dnorm(beta.prop[k], mean = 0, sd = sd.beta[k], log = TRUE)
    # log-proposal density at beta.ini[k] and beta.prop[k]
    logprop.ini                  <- dnorm(beta.ini[k],  mean = beta.prop.mean,     sd = sqrt(beta.prop.variance),     log = TRUE)
    logprop.prop                 <- dnorm(beta.prop[k], mean = beta.prop.mean.ini, sd = sqrt(beta.prop.variance.ini), log = TRUE)
    
    logR                         <- loglik.prop - loglik.ini + logprior.prop - logprior.ini + logprop.ini - logprop.prop
    
    u = log(runif(1)) < logR
    
    if(u == 1){
      beta.ini[k] <- beta.prop[k]
      Xbeta	      <- Xbeta.prop
    }
    
    accept[k]     <- accept[k] + u
    
  } 
  
  list(beta.ini = beta.ini, accept = accept, Xbeta = Xbeta)
  
}
#'
UpdateBH.lambda.SMCFM       <- function(survObj, priorPara, mcmcPara, ini){
  
  X            <- survObj$X
  Z            <- survObj$Z
  y            <- survObj$t
  delta        <- survObj$di
  J            <- priorPara$J
  v            <- priorPara$v
  s            <- priorPara$s
  b.ini        <- ini$b.ini
  beta.ini     <- ini$beta.ini
  a            <- priorPara$a 
  b            <- priorPara$b 
  lambda.ini   <- ini$lambda.ini
  theta.ini    <- ini$theta.ini
  n            <- nrow(X)
  p            <- length(lambda.ini)
  accept       <- rep(0, p)
  
  # Compute cumulative hazard and hazard of each t[i]
  H0_h0                    <- cumulative_hazard(y, interval_bounds=s, lambda=lambda.ini)
  H0                       <- H0_h0$H0
  h0                       <- H0_h0$h0 
  
  Xbeta                    <- X %*% beta.ini
  exp.Xbeta                <- exp(Xbeta)
  Zb                       <- Z %*% b.ini 
  exp.Zb                   <- exp(Zb)
  
  for (k in 1:J) {
    H1                     <- H0 * exp.Xbeta 
    H20                    <- H1 / theta.ini  
    H2                     <- 1 + H20         
    term1                  <- ( log(h0) - (theta.ini+1) * log1p(H20) ) * delta 
    term2                  <- log1p( exp.Zb * H2^(-theta.ini) ) * (1-delta) 
    loglik.ini             <- sum(term1 + term2)
    # proposal value
    lambda.prop.shape      <- 1 *lambda.ini[k]  
    lambda.prop.rate       <- 1
    lambda.prop            <- lambda.ini
    lambda.prop[k]	       <- mean( rgamma(5, shape = lambda.prop.shape, rate =  lambda.prop.rate ) )
    while(lambda.prop[k] < 0.0001){
      lambda.prop[k]	     <- mean( rgamma(5, shape = lambda.prop.shape, rate =  lambda.prop.rate ) )
    }
    # Calculating acceptance probability 
    
    # baseline cumulative hazard function H0 updated lambda[k]!
    H0_h0.prop                <- cumulative_hazard(y, interval_bounds = s, lambda = lambda.prop)
    H0.prop                   <- H0_h0.prop$H0 
    h0.prop                   <- H0_h0.prop$h0 
    H1.prop                   <- H0.prop * exp.Xbeta 
    H20.prop                  <- H1.prop / theta.ini  
    H2.prop                   <- 1 + H20.prop  
    term1.prop                <- ( log(h0.prop) - (theta.ini+1)* log1p(H20.prop) ) * delta 
    term2.prop                <- log1p( exp.Zb * H2.prop^(-theta.ini) ) * (1-delta)
    loglik.prop               <- sum(term1.prop + term2.prop)                   # log-likelihood for lambda_k at initial value
    lambda.prop.shape.ini     <- 1*lambda.prop[k]  
    lambda.prop.rate.ini      <- 1 
    
    # log-prior density at lambda.ini[k] and lambda.prop[k]
    logprior.ini              <- dgamma(lambda.ini[k],  shape = a , rate = b, log = TRUE)  
    logprior.prop             <- dgamma(lambda.prop[k], shape = a , rate = b, log = TRUE) 
    # log-proposal density at lambda.ini[k] and lambda.prop[k] from the proposal density 
    logprop.ini               <- dgamma(lambda.ini[k],  shape = lambda.prop.shape,     rate = lambda.prop.rate,     log = TRUE)
    logprop.prop              <- dgamma(lambda.prop[k], shape = lambda.prop.shape.ini, rate = lambda.prop.rate.ini, log = TRUE)
    
    logR                      <- loglik.prop - loglik.ini + logprior.prop - logprior.ini + logprop.ini - logprop.prop
    
    u = log(runif(1)) < logR
    
    if(u == 1){
      lambda.ini[k]          <- lambda.prop[k]              
      H0	                   <- H0.prop
      h0                     <- h0.prop
    }
    
    accept[k]                <- accept[k] + u
  } 
  
  list(lambda.ini = lambda.ini, accept = accept, H0 = H0)
  
}
#'
Update.theta.SMCFM          <- function(survObj, priorPara, mcmcPara, ini){
  
  X                <- survObj$X
  Z                <- survObj$Z
  y                <- survObj$t
  delta            <- survObj$di
  J                <- priorPara$J
  v                <- priorPara$v
  s                <- priorPara$s
  b.ini            <- ini$b.ini
  beta.ini         <- ini$beta.ini
  c                <- priorPara$c 
  d                <- priorPara$d 
  lambda.ini       <- ini$lambda.ini
  theta.ini        <- ini$theta.ini
  n                <- nrow(X)
  accept           <- rep(0, length(theta.ini)) 
  
  H0_h0            <- cumulative_hazard(y, interval_bounds=s, lambda=lambda.ini)
  H0               <- H0_h0$H0 
  h0               <- H0_h0$h0 
  Xbeta            <- X %*% beta.ini
  exp.Xbeta        <- exp(Xbeta)
  Zb               <- Z %*% b.ini
  exp.Zb           <- exp(Zb)
  H1               <- H0 * exp.Xbeta
  H20              <- H1 / theta.ini 
  H2               <- 1 + H20
  term1            <- -(theta.ini+1)* log1p(H20) * delta 
  term2            <- log1p(exp.Zb * H2^(-theta.ini) ) * (1-delta) 
  loglik.ini       <- sum(term1 + term2) 
  # proposal value
  theta.prop.shape <- 1 *theta.ini
  theta.prop.rate  <- 1 
  
  theta.prop       <- theta.ini
  theta.prop    	 <- mean( rgamma(5, shape = theta.prop.shape, rate =  theta.prop.rate ) ) 
  while(theta.prop < 0.0001){
    theta.prop	   <- mean( rgamma(5, shape = theta.prop.shape, rate =  theta.prop.rate ) )
  }
  # Calculating acceptance probability 
  H20.prop         <- H1 / theta.prop 
  H2.prop          <- 1 + H20.prop    
  term1.prop       <- -(theta.prop+1)* log1p(H20.prop) * delta 
  term2.prop       <- log1p(exp.Zb * H2.prop^(-theta.ini) ) * (1-delta)  
  loglik.prop      <- sum(term1.prop + term2.prop)  
  
  theta.prop.shape.ini<- 1*theta.prop 
  theta.prop.rate.ini <- 1 
  # log-prior density at theta.ini and theta.prop
  logprior.ini        <- dgamma(theta.ini,  shape = c , rate = d, log = TRUE)  
  logprior.prop       <- dgamma(theta.prop, shape = c , rate = d, log = TRUE) 
  # log-proposal density at theta.ini and theta.prop from the proposal density which should be positive!
  logprop.ini         <- dgamma(theta.ini,  shape = theta.prop.shape,     rate = theta.prop.rate,     log = TRUE)
  logprop.prop        <- dgamma(theta.prop, shape = theta.prop.shape.ini, rate = theta.prop.rate.ini, log = TRUE)
  
  logR                <- loglik.prop - loglik.ini + logprior.prop - logprior.ini + logprop.ini - logprop.prop
  
  u = log(runif(1)) < logR
  if(u == 1){
    theta.ini       <- theta.prop 
  }
  accept            <- accept + u
  
  list(theta.ini = theta.ini, accept = accept )
  
}