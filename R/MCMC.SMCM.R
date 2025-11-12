#' @title Function to MCMC algorithm for Semi-parametric Mixture Cure Model (SMCM)
#'
#' @description
#' This function obtains posterior samples for all model parameters using an MCMC algorithm 
#' for a single chain when applying the \code{MCMC.SMCM} framework in Kızılaslan and Vitelli (2025) to the data.
#'
#' @name MCMC.SMCM
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar flush.console
#'
#' @param survObj a list containing observed data from \code{n} subjects,
#' \code{t}, \code{di}, \code{X}, \code{Z}. See details for more information.
#' @param priorPara a list containing prior parameter values for the SMCM model; among c("r1", "delta1", "r2", "delta2", "a", "b", "s", "J").
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
#'    \item{mcmcOutcome}{The list containing posterior samples for the remaining model parameters.}
#'    \item{t}{Processed times to the event data.}
#'    \item{di}{Processed censoring indicators for the event time.}
#'    \item{loglik.matrix}{Matrix of the log-likelihood of the model over MCMC samples when \code{save_loglik = TRUE}.}
#' 
#' @references 
#' Kızılaslan, F., Vitelli, V. (2025). Bayesian Semiparametric Mixture Cure (Frailty) Model, arXiv \url{https://arxiv.org/???}.
#' 
#' @export
MCMC.SMCM                   <- function(survObj, priorPara, initial, nIter = 6000, save_loglik = 0){
  
  r1		       	            <- priorPara$r1
  delta1			              <- priorPara$delta1
  r2		       	            <- priorPara$r2
  delta2			              <- priorPara$delta2
  a                         <- priorPara$a
  b                         <- priorPara$b
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
  tauSq.p						        <- tauSq
  tauStarSq.p						    <- tauStarSq
  mcmcOutcome$sigmaSq.p	    <- sigmaSq
  mcmcOutcome$sigmaStarSq.p	<- sigmaStarSq
  mcmcOutcome$etaSq.p		    <- etaSq
  mcmcOutcome$etaStarSq.p	  <- etaStarSq
  mcmcOutcome$accept.beta	  <- c( rep(0, length(beta.ini)) )
  mcmcOutcome$accept.b		  <- c( rep(0, length(b.ini)) )
  mcmcOutcome$accept.lambda	<- c( rep(0, length(lambda.ini)) )
  
  if (save_loglik == 1) {
    loglik.p                  <- logLik_SMCM(survObj, priorPara, b.coef = b.p, beta.coef = beta.p, lambda.coef = lambda.p ) }     # at initial values
  
  pb  <- txtProgressBar(min = 0, max = nIter, style = 3)
  
  # MCMC sampling
  for(M in 1:nIter){
   
    # Updating regression parameters b_k
    sampleRP.b	                    <- UpdateRP.b(survObj, priorPara, mcmcPara, ini)
    b.ini	    <- ini$b.ini	        <- sampleRP.b$b.ini
    Zb		    <- ini$Zb	            <- sampleRP.b$Zb
    mcmcOutcome$accept.b	          <- mcmcOutcome$accept.b + sampleRP.b$accept
    
    # Updating regression parameters beta_k:
    sampleRP.beta	                  <- UpdateRP.beta(survObj, priorPara, mcmcPara, ini)
    beta.ini	 <- ini$beta.ini	    <- sampleRP.beta$beta.ini
    Xbeta		   <- ini$Xbeta	        <- sampleRP.beta$Xbeta
    mcmcOutcome$accept.beta	        <- mcmcOutcome$accept.beta + sampleRP.beta$accept
    
    # Updating baseline hazard paramaters lambda_k: 
    sampleBH.lambda	                <- UpdateBH.lambda(survObj, priorPara, mcmcPara, ini)
    lambda.ini	                    <- ini$lambda.ini	<- sampleBH.lambda$lambda.ini
    mcmcOutcome$accept.lambda	      <- mcmcOutcome$accept.lambda + sampleBH.lambda$accept
    
    # Updating 1/tauSq
    tauSq	     <- ini$tauSq	        <- UpdateTau(survObj, priorPara, mcmcPara, ini)
    
    # Updating 1/tauStarSq
    tauStarSq	 <- ini$tauStarSq	    <- UpdateTauStar(survObj, priorPara, mcmcPara, ini)
    
    # Updating sigmaSq
    sigmaSq	   <- ini$sigmaSq	      <- UpdateSigma(survObj, priorPara, mcmcPara, ini)
    
    # Updating sigmaStarSq
    sigmaStarSq	<- ini$sigmaStarSq  <- UpdateSigmaStar(survObj, priorPara, mcmcPara, ini)
    
    # Updating etaSq		
    etaSq	      <- ini$etaSq	      <- UpdateEta(survObj, priorPara, mcmcPara, ini)
    
    # Updating etaStarSq		
    etaStarSq	  <- ini$etaStarSq	  <- UpdateEtaStar(survObj, priorPara, mcmcPara, ini)
    
    ini$sd.beta <- sd.beta	        <- sqrt(sigmaStarSq*tauStarSq) 
    ini$sd.b	  <- sd.b             <- sqrt(sigmaSq*tauSq)
    
    if (save_loglik == 1) {
      loglik.ini                    <- logLik_SMCM(survObj, priorPara, b.coef = b.ini, beta.coef = beta.ini, lambda.coef = lambda.ini ) 
    }
    
    ###### storing posterior samples
    beta.p						              <- rbind(beta.p, beta.ini, deparse.level = 0)
    b.p						                  <- rbind(b.p, b.ini, deparse.level = 0)
    lambda.p				                <- rbind(lambda.p, lambda.ini, deparse.level = 0)
    if (save_loglik == 1) { 
      loglik.p                      <- rbind(loglik.p, loglik.ini, deparse.level = 0) 
    } 
    
    tauSq.p						             <- rbind(tauSq.p, tauSq, deparse.level = 0)
    tauStarSq.p						         <- rbind(tauStarSq.p, tauStarSq, deparse.level = 0)
    mcmcOutcome$sigmaSq.p	         <- c(mcmcOutcome$sigmaSq.p, sigmaSq)
    mcmcOutcome$sigmaStarSq.p	     <- c(mcmcOutcome$sigmaStarSq.p, sigmaStarSq)
    mcmcOutcome$etaSq.p		         <- c(mcmcOutcome$etaSq.p, etaSq)
    mcmcOutcome$etaStarSq.p	       <- c(mcmcOutcome$etaStarSq.p, etaStarSq)
    mcmcOutcome$ini				         <- ini
    
    # update progress bar every 10% of iterations
    if (M %% (nIter / 10) == 0 || M == nIter) {
      setTxtProgressBar(pb, M)
      flush.console()
    }
    
  } 
  
  close(pb)
  
  if (save_loglik == 1) {
    ret      <- list(beta.p = beta.p, b.p = b.p, tauSq.p = tauSq.p,  tauStarSq.p = tauStarSq.p, lambda.p = lambda.p,  mcmcOutcome = mcmcOutcome, t=survObj$t, di=survObj$di, loglik.matrix = loglik.p ) 
  }else{
    ret      <- list(beta.p = beta.p, b.p = b.p, tauSq.p = tauSq.p,  tauStarSq.p = tauStarSq.p, lambda.p = lambda.p,  mcmcOutcome = mcmcOutcome, t=survObj$t, di=survObj$di  )
  }
  return(ret)
} 
#'
UpdateRP.b  <- function(survObj, priorPara, mcmcPara, ini){
  
  Z         <- survObj$Z
  X         <- survObj$X
  y         <- survObj$t
  delta     <- survObj$di
  J         <- priorPara$J
  v         <- priorPara$v
  s         <- priorPara$s
  b.ini     <- ini$b.ini	 
  beta.ini  <- ini$beta.ini
  lambda.ini<- ini$lambda.ini
  sigmaSq   <- ini$sigmaSq	
  tauSq     <- ini$tauSq	
  sd.b      <- ini$sd.b	
  n         <- nrow(Z)
  p         <- length(b.ini)
  accept    <- rep(0, p)
  
  # Compute cumulative hazard and hazard for each t[i]
  H0_h0     <- cumulative_hazard(y, interval_bounds = s, lambda = lambda.ini)
  H0        <- H0_h0$H0 
  h0        <- H0_h0$h0 
  Zb        <- Z %*% b.ini
  Xbeta     <- X %*% beta.ini
  exp.Xbeta <- exp(Xbeta)
  H1        <- H0 * exp.Xbeta 
  
  for (k in 1:p) {
    
    result  <- NULL
    while ( is.null(result) ) {
      result <- tryCatch({
        
        exp.Zb                  <- exp(Zb)
        term1                   <- (Zb * delta) - log1p(exp.Zb)
        term2                   <- log1p( exp(-H1 + Zb) ) * (1-delta)  
        loglik.ini              <- sum(term1 + term2)                           # log-likelihood for b_k at initial value
        # proposal value
        b.prop.mean             <- b.ini[k] 
        b.prop.variance         <- 1^2 
        b.prop                  <- b.ini
        b.prop[k]	              <- mean( rnorm(5, mean = b.prop.mean, sd = sqrt(b.prop.variance) ) )
        # Calculating acceptance probability
        Zb.prop	                <- Zb - Z[,k] * b.ini[k] + Z[,k] * b.prop[k]    # updating Zb with new b[k]
        exp.Zb.prop	            <- exp(Zb.prop)
        term1.prop              <- (Zb.prop * delta) - log1p(exp.Zb.prop)
        term2.prop              <- log1p( exp(-H1 + Zb.prop) ) * (1-delta)  
        loglik.prop             <- sum(term1.prop + term2.prop )                # log-likelihood at proposal value of b[k]
        b.prop.mean.ini         <- b.prop[k] 
        b.prop.variance.ini     <- 1^2
        #log-prior density at b.ini[k] and b.prop[k]
        logprior.ini            <- dnorm(b.ini[k],  mean = 0 , sd = sd.b[k], log = TRUE) 
        logprior.prop           <- dnorm(b.prop[k], mean = 0 , sd = sd.b[k], log = TRUE)
        # log-proposal density at b.ini[k] and b.prop[k]
        logprop.ini             <- dnorm(b.ini[k], mean = b.prop.mean, sd = sqrt(b.prop.variance), log = TRUE)
        logprop.prop            <- dnorm(b.prop[k], mean = b.prop.mean.ini, sd = sqrt(b.prop.variance.ini), log = TRUE)
        
        logR                    <- loglik.prop - loglik.ini + logprior.prop - logprior.ini + logprop.ini - logprop.prop
        if ( is.nan(logR) ) stop("NaN encountered!")  
        logR
      }, 
      error = function(e) { 
        cat("Error occurred in UpdateRP.b, retrying...\n")
        NULL  
      })
    } 
    
    u = log(runif(1)) < logR
    if(u == 1){
      b.ini[k]        <- b.prop[k]
      Zb	            <- Zb.prop
    }
    
    accept[k]        <- accept[k] + u
    
  } 
  
  list(b.ini = b.ini, accept = accept, Zb = Zb)
}
#'
UpdateRP.beta <- function(survObj, priorPara, mcmcPara, ini){
  
  X           <- survObj$X
  Z           <- survObj$Z
  y           <- survObj$t
  delta       <- survObj$di
  J           <- priorPara$J
  v           <- priorPara$v
  s           <- priorPara$s
  b.ini       <- ini$b.ini
  beta.ini    <- ini$beta.ini	 
  sigmaStarSq <- ini$sigmaStarSq	
  tauStarSq   <- ini$tauStarSq
  lambda.ini  <- ini$lambda.ini
  sd.beta     <- ini$sd.beta	
  n           <- nrow(X)
  p           <- length(beta.ini)
  accept      <- rep(0, p)
  Xbeta       <- X %*% beta.ini 
  Zb          <- Z %*% b.ini 
  
  H0_h0    <- cumulative_hazard(y, interval_bounds = s, lambda = lambda.ini)
  H0       <- H0_h0$H0 
  h0       <- H0_h0$h0 
  
  for (k in 1:p) {
    
    result <- NULL
    while ( is.null(result) ) {
      result <- tryCatch({
        
        exp.Xbeta                   <- exp(Xbeta)
        H1                          <- H0 * exp.Xbeta 
        term1                       <- (Xbeta - H1) * delta 
        term2                       <- log1p(exp(-H1 + Zb) )* (1-delta)
        loglik.ini                  <- sum(term1 + term2)                       # log-likelihood for beta_k at initial value
        beta.prop.mean              <- beta.ini[k] 
        beta.prop.variance          <- 1^2 
        beta.prop                   <- beta.ini
        beta.prop[k]	              <- mean( rnorm(5, mean = beta.prop.mean, sd = sqrt(beta.prop.variance) ) )
        # Calculating acceptance probability
        Xbeta.prop	                 <- Xbeta - X[,k] * beta.ini[k] + X[,k] * beta.prop[k] # updating Xbeta with new beta[k]
        exp.Xbeta.prop	             <- exp(Xbeta.prop)
        H1.prop                      <- H0 * exp.Xbeta.prop 
        term1.prop                   <- (Xbeta.prop - H1.prop) * delta 
        term2.prop                   <- log1p( exp(-H1.prop + Zb) ) * (1-delta)  
        loglik.prop                  <- sum(term1.prop + term2.prop)            # log-likelihood  at proposal value of beta[k]
        beta.prop.mean.ini           <- beta.prop[k] 
        beta.prop.variance.ini       <- 1^2 
        #log-prior density at beta.ini[k] and beta.prop[k]
        logprior.ini                 <- dnorm(beta.ini[k],  mean = 0, sd = sd.beta[k], log = TRUE) 
        logprior.prop                <- dnorm(beta.prop[k], mean = 0, sd = sd.beta[k], log = TRUE)
        # log-proposal density at beta.ini[k] and beta.prop[k]
        logprop.ini                  <- dnorm(beta.ini[k],  mean = beta.prop.mean,     sd = sqrt(beta.prop.variance),     log = TRUE)
        logprop.prop                 <- dnorm(beta.prop[k], mean = beta.prop.mean.ini, sd = sqrt(beta.prop.variance.ini), log = TRUE)
        
        logR                         <- loglik.prop - loglik.ini + logprior.prop - logprior.ini + logprop.ini - logprop.prop
        if ( is.nan(logR) ) stop("NaN encountered!")  
        logR
      }, 
      error = function(e) { 
        cat("Error occurred in UpdateRP.beta, retrying...\n")
        NULL  # Return NULL to trigger retry
      })
    } 
    
    u = log(runif(1)) < logR
    
    if(u == 1){
      beta.ini[k]    <- beta.prop[k]
      Xbeta	         <- Xbeta.prop
    }
    
    accept[k]        <- accept[k] + u
    
  } 
  
  list(beta.ini = beta.ini, accept = accept, Xbeta = Xbeta)
  
}
#'
UpdateBH.lambda <- function(survObj, priorPara, mcmcPara, ini){
  X             <- survObj$X
  Z             <- survObj$Z
  y             <- survObj$t
  delta         <- survObj$di
  J             <- priorPara$J
  v             <- priorPara$v
  s             <- priorPara$s
  b.ini         <- ini$b.ini
  beta.ini      <- ini$beta.ini
  a             <- priorPara$a 
  b             <- priorPara$b 
  lambda.ini    <- ini$lambda.ini
  n             <- nrow(X)
  p             <- length(lambda.ini)
  accept        <- rep(0, p) 
  
  H0_h0         <- cumulative_hazard(y, interval_bounds = s, lambda = lambda.ini)
  H0            <- H0_h0$H0 
  h0            <- H0_h0$h0 
  Xbeta         <- X %*% beta.ini
  exp.Xbeta     <- exp(Xbeta)
  Zb            <- Z %*% b.ini
  
  for (k in 1:p) {
    
    result <- NULL
    while ( is.null(result) ) {
      result <- tryCatch({
        
        H1                    <- H0 * exp.Xbeta
        term1                 <- (log(h0) - H1) * delta
        term2                 <- log1p(exp(-H1 + Zb) )* (1-delta)
        loglik.ini            <- sum(term1 + term2)                             # log-likelihood for lambda_k at initial value
        # proposal value
        lambda.prop.shape     <- 1 *lambda.ini[k]  
        lambda.prop.rate      <- 1                                              # variance = shape/ rate^2
        
        lambda.prop           <- lambda.ini
        lambda.prop[k]	      <- mean( rgamma(5, shape = lambda.prop.shape, rate =  lambda.prop.rate ) ) 
        while(lambda.prop[k]  < 0.0001){
          lambda.prop[k]	    <- mean( rgamma(5, shape = lambda.prop.shape, rate =  lambda.prop.rate ) )
        }
        # baseline cumulative hazard function H0 for updated lambda[k]!
        H0_h0.prop            <- cumulative_hazard(y, interval_bounds = s, lambda = lambda.prop)
        H0.prop               <- H0_h0.prop$H0 
        h0.prop               <- H0_h0.prop$h0 
        H1.prop               <- H0.prop * exp.Xbeta
        term1.prop            <- (log(h0.prop) - H1.prop) * delta
        term2.prop            <- log1p(exp(-H1.prop + Zb) )* (1-delta)
        loglik.prop           <- sum(term1.prop + term2.prop)                   # log-likelihood for lambda_k at initial value
        lambda.prop.shape.ini <- 1*lambda.prop[k]  
        lambda.prop.rate.ini  <- 1                                              # variance=shape/ rate^2
        # log-prior density at lambda.ini[k] and lambda.prop[k]
        logprior.ini          <- dgamma(lambda.ini[k],  shape = a , rate = b, log = TRUE)  
        logprior.prop         <- dgamma(lambda.prop[k], shape = a , rate = b, log = TRUE) 
        # log-proposal density at lambda.ini[k] and lambda.prop[k] from the proposal density 
        logprop.ini           <- dgamma(lambda.ini[k],  shape = lambda.prop.shape,     rate = lambda.prop.rate,     log = TRUE)
        logprop.prop          <- dgamma(lambda.prop[k], shape = lambda.prop.shape.ini, rate = lambda.prop.rate.ini, log = TRUE)
        logR                  <- loglik.prop - loglik.ini + logprior.prop - logprior.ini + logprop.ini - logprop.prop
        
        if ( is.nan(logR) ) stop("NaN encountered!")  
        logR
      }, 
      error = function(e) { 
        cat("Error occurred in UpdateBH.lambda, retrying...\n")
        NULL  
      })
    } 
    u = log(runif(1))        < logR
    
    if(u == 1){
      lambda.ini[k]           <- lambda.prop[k] #updated kth lambda and H0!
      H0	                    <- H0.prop
      h0                      <- h0.prop
    }
    accept[k]                 <- accept[k] + u
  } 
  
  list(lambda.ini = lambda.ini, accept = accept, H0 = H0)
  
}