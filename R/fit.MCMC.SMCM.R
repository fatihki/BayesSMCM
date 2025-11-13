#' fit.MCMC.SMCM
#' @title Fit Bayesian Semi-parametric Mixture Cure Model (SMCM) via MCMC algorithm
#'
#' @description
#' This is the main function to fit the Bayesian SMCM model using MCMC algorithm in Kızılaslan and Vitelli (2025).
#'
#' @name fit.MCMC.SMCM
#' 
#' @import parallel
#' @import doParallel
#' @importFrom foreach foreach  %dopar% 
#'
#' @param data a list containing observed data from \code{n} subjects,
#' \code{observed_time}, \code{delta}, \code{X}, \code{Z}. See details for more information.
#' @param hyperpar a list containing prior parameter values for the SMCM model; among c("r1", "delta1", "r2", "delta2", "a", "b").
#' See details for more information.
#' @param nchains Number of MCMC chains.
#' @param nIter Number of MCMC iterations.
#' @param warmup Number of warmup iterations.
#' @param thin Thinning interval for MCMC samples.
#' @param mcmc.parallel Parallelization method to run nchains for MCMC algorithm, can be \code{lapply}, 
#' \code{foreach} from \pkg{doParallel} package or \code{parLapply} from \pkg{parallel} package.
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
#' @return A list containing:
#' \item{fit}{A list containing all results across MCMC chains.}
#' \item{fit.results}{Summary of posterior estimation results, including posterior means, credible intervals, and 
#' convergence diagnostics for all model parameters.}
#' \item{b.chains}{Posterior samples for the regression coefficients \eqn{\bm b} across all MCMC chains.}
#' \item{beta.chains}{Posterior samples for the regression coefficients \eqn{\bm \beta} across all MCMC chains.}
#' \item{lambda.chains}{Posterior samples for the regression coefficients \eqn{\bm \lambda} across all MCMC chains.}
#' \item{Rhat}{Gelman–Rubin convergence statistic for assessing MCMC chain convergence.}
#' \item{loo.results}{Model fit statistics including DIC and LPML when \code{save_loglik} = TRUE.}
#' \item{loglik.matrix.chains}{Matrix of the log-likelihood of the model over MCMC samples when \code{save_loglik = TRUE}.}
#' \item{data}{Processed survival data used for model fitting.}
#' \item{priors}{Prior parameter settings used in the model.}
#' \item{mcmc.info}{nchains, nIter, warmup and thin are used in the model.}
#' 
#' @examples
#' # Fit SMCM model to data generated from Scenario 1 in our simulation study
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
#' nchains = 2; nIter = 2500; warmup = 500; thin = 10
#' priorPar = list( r1 = 1, delta1 = 0.0001, r2 = 1, delta2 = 0.0001, a = 0.1, b = 0.1 )
#' out.smcm.mcmc = fit.MCMC.SMCM ( data = dat1, hyperpar = priorPar, nchains, nIter, warmup, 
#'   thin, mcmc.parallel = "parLapply", standardize = FALSE, probs = 0, save_loglik = 1, seed = 2025 )
#' print_smcm( out.smcm.mcmc, stan.model = FALSE, frailty = FALSE )
#' 
#' @references 
#' Kızılaslan, F., Vitelli, V. (2025). Bayesian Semiparametric Mixture Cure (Frailty) Model, arXiv \url{https://arxiv.org/???}.
#' 
#' @export
fit.MCMC.SMCM            <- function(data, hyperpar = list(), nchains = 2, nIter = 6000, warmup = 1000, thin = 20, mcmc.parallel =  "parLapply",
                                     standardize,  probs, save_loglik, seed = NULL  ){
  
  if (!is.null(seed)) set.seed(seed)

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

  initial               <- list()
  initial$beta.ini      <- rep(1, ncol(data$X)  )                      
  initial$b.ini         <- rep(1, ncol(data$Z)  )                         
  initial$sigmaSq       <- runif(1, 0.1, 10) 
  initial$sigmaStarSq   <- runif(1, 0.1, 10) 
  initial$etaSq         <- 1  
  initial$etaStarSq     <- 1  
  initial$tauSq         <- rexp( length(initial$b.ini), rate = initial$etaSq/2) 
  initial$tauStarSq     <- rexp( length(initial$beta.ini), rate = initial$etaStarSq/2) 
  initial$lambda.ini    <- rgamma( priorPara$J, 1, 1)
    
    if( nchains > 1){
      
      if (!is.null(seed))  seeds <- sample(1:1e6, nchains)                      # random seed generated for per chain if we use fixed seed
      
      if(mcmc.parallel == "lapply"){
        fun <- function(i) {
          if (!is.null(seed)) set.seed(seeds[i])
          # update initial points for each chain
          initial$beta.ini <- rep(0.5*i, ncol(survObj$X) )
          initial$b.ini    <- rep(0.5*i, ncol(survObj$Z) )
          return( MCMC.SMCM(survObj, priorPara, initial,  nIter, save_loglik ) )
        }
        start.mcmc.smc.time   <- Sys.time()
        fit                   <- lapply( 1:nchains, FUN=fun)
        mcmc.smc.time         <- as.numeric(Sys.time() - start.mcmc.smc.time, units = "mins")
      }
      
      if(mcmc.parallel == "foreach"){
        start.mcmc.smc.time   <- Sys.time()
        cl1                   <- parallel::makeCluster( nchains )
        doParallel::registerDoParallel(cl1)
        # parallel::clusterCall(cl1, function() library(doParallel) )
        # parallel::clusterCall(cl1, function() library(BayesSMCM) )
        parallel::clusterExport(cl1, varlist = "MCMC.SMCM", envir = asNamespace("BayesSMCM"))
        fit                <- foreach( i=1:nchains) %dopar% {
          if (!is.null(seed)) set.seed(seeds[i])
          initial$beta.ini <- rep(0.5*i, ncol(survObj$X) )
          initial$b.ini    <- rep(0.5*i, ncol(survObj$Z) )
          fit.mcmc         <- MCMC.SMCM(survObj, priorPara, initial,  nIter,  save_loglik )
          return(fit.mcmc)
        }
        on.exit(stopCluster(cl1))
        mcmc.smc.time      <- as.numeric(Sys.time() - start.mcmc.smc.time, units = "mins")
      }
      
      if(mcmc.parallel == "parLapply"){
        fun <- function(i, survObj, priorPara, initial, nIter, save_loglik  ) {
          if (!is.null(seed)) set.seed(seeds[i])
          initial$beta.ini <- rep(0.5*i, ncol(survObj$X) )
          initial$b.ini    <- rep(0.5*i, ncol(survObj$Z) )
          return( MCMC.SMCM(survObj, priorPara, initial,  nIter,  save_loglik ) )
        }
        start.mcmc.smc.time <- Sys.time()
        cl1                 <- parallel::makeCluster( nchains )
        # parallel::clusterCall(cl1, function() library(doParallel) )
        # parallel::clusterCall(cl1, function() library(parallel) )
        # parallel::clusterCall(cl1, function() library(BayesSMCM) )
        parallel::clusterExport(cl1, varlist = "MCMC.SMCM", envir = asNamespace("BayesSMCM"))
        fit                  <- parallel::parLapply(cl1, 1:nchains, fun, survObj, priorPara, initial, nIter,  save_loglik )
        on.exit(parallel::stopCluster(cl1))
        mcmc.smc.time       <- as.numeric(Sys.time() - start.mcmc.smc.time, units = "mins")
      }
      
      Rhat                  <- Rhat.mcmc(fit, warmup, thin )                   
      b.mcmc.chains         <- mcmc.burn.thin(fit, warmup, thin, par = "b")     # mcmc chains after removed warmup iterations and thinning
      b.hat.mcmc            <- colMeans( matrix( unlist( lapply(b.mcmc.chains, function(x) colMeans(x) ) ), nrow = nchains, byrow = T) )
      
      beta.mcmc.chains      <- mcmc.burn.thin(fit, warmup, thin, par = "beta")
      if( ncol( as.matrix(beta.mcmc.chains[[1]]) ) == 1){
        beta.hat.mcmc       <- colMeans( matrix( unlist( lapply(beta.mcmc.chains, function(x) mean(x) ) ),  nrow = nchains, byrow = T) )
      }else{
        beta.hat.mcmc       <- colMeans( matrix( unlist( lapply(beta.mcmc.chains, function(x) colMeans(x) ) ),  nrow = nchains, byrow = T) )
      }
      
      lambda.mcmc.chains    <- mcmc.burn.thin(fit, warmup, thin, par = "lambda")
      if( ncol( as.matrix(lambda.mcmc.chains[[1]]) ) == 1){
        lambda.hat.mcmc     <- colMeans( matrix( unlist( lapply(lambda.mcmc.chains, function(x) mean(x) ) ),   nrow = nchains, byrow = T) ) 
      }else{
        lambda.hat.mcmc     <- colMeans( matrix( unlist( lapply(lambda.mcmc.chains, function(x) colMeans(x) ) ), nrow = nchains, byrow = T) ) 
      }
      if(save_loglik == 1){
        loglik.matrix.chains <- mcmc.burn.thin(fit, warmup, thin, par = "loglik")
        loo.results          <- mcmc.dic.lpml( survObj, priorPara, loglik.matrix.chains, b.hat.mcmc, beta.hat.mcmc, lambda.hat.mcmc, stan.model = FALSE, frailty = FALSE )
      }
    }else{
      # nchains = 1 case
      start.mcmc.smc.time   <- Sys.time()
      fit                   <- MCMC.SMCM( survObj, priorPara, initial, nIter,  save_loglik )
      mcmc.smc.time         <- as.numeric(Sys.time() - start.mcmc.smc.time, units = "mins")
      
      # apply warmup and thin to all mcmc iterations
      b.mcmc.chains         <- fit$b.p[ -c(1:(warmup+1)), ][seq(1, nIter - warmup, by = thin), ] #even if we have only one chain, we use "b.mcmc.chains" to keep the variable name same
      b.hat.mcmc            <- colMeans( b.mcmc.chains ) # colMeans( fit$b.p[-c(1:(warmup+1)), ][seq(1, nIter - warmup, by = thin), ]  )
      
      if( ncol( as.matrix(fit$beta.p) ) == 1){
        beta.mcmc.chains    <- fit$beta.p[ -c(1:(warmup+1)) ][seq(1, nIter - warmup, by = thin) ]
        beta.hat.mcmc       <- mean( beta.mcmc.chains ) # mean( fit$beta.p[-c(1:(warmup+1)) ][seq(1, nIter - warmup, by = thin) ]  )
      }else{
        beta.mcmc.chains    <- fit$beta.p[ -c(1:(warmup+1)), ][seq(1, nIter - warmup, by = thin), ]  
        beta.hat.mcmc       <- colMeans( beta.mcmc.chains ) # colMeans( fit$beta.p[-c(1:(warmup+1)), ][seq(1, nIter - warmup, by = thin), ]  ) }
      }
      
      if( ncol( as.matrix(fit$lambda.p) ) == 1){
        lambda.mcmc.chains  <- fit$lambda.p[ -c(1:(warmup+1)) ][ seq(1, nIter - warmup, by = thin) ]
        lambda.hat.mcmc     <- mean( lambda.mcmc.chains ) # mean( fit$lambda.p[-c(1:(warmup+1)) ][seq(1, nIter - warmup, by = thin) ]  )
      }else{
        lambda.mcmc.chains  <- fit$lambda.p[-c(1:(warmup+1)), ][seq(1, nIter - warmup, by = thin), ] 
        lambda.hat.mcmc     <- colMeans( lambda.mcmc.chains ) # colMeans( fit$lambda.p[-c(1:(warmup+1)), ][seq(1, nIter - warmup, by = thin), ]  )
      }
      if(save_loglik == 1){
        loglik.matrix.chains <- fit$loglik.matrix[ -c(1:(warmup+1)), ][ seq(1, nIter - warmup, by = thin), ] 
        loo.results          <- mcmc.dic.lpml( survObj, priorPara, loglik.matrix.chains, b.hat.mcmc, beta.hat.mcmc, lambda.hat.mcmc, stan.model = FALSE, frailty = FALSE )
      }
      Rhat                  <- NULL 
    }
    
    pz.hat.mcmc              <- logit( as.matrix(survObj$Z), b.hat.mcmc) 
    fit.results              <- c( b.hat = b.hat.mcmc , beta.hat =  beta.hat.mcmc, 
                                   lambda.hat =  lambda.hat.mcmc,  pz.hat = mean(pz.hat.mcmc),
                                   time = mcmc.smc.time )
 
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
    
    class(ret) <- "SMCM"
    
    return(ret)   
}
