#' Internal helper functions for MCMC algorithm for both SMCM and SMCFM 
#'
#' @importFrom SuppDists rinvGauss
#' @importFrom LearnBayes rigamma
#' 
#' @noRd
cumulative_hazard <- function(t, interval_bounds, lambda) {
  # t: vector of observed times
  # interval_bounds (s): vector of interval breakpoints (length J+1) 
  # lambda: vector of hazard rates for each interval (length J)
  
  n                 <- length(t)              # Number of subjects
  J                 <- length(lambda)         # Number of intervals
  
  cumulative_hazard <- c()
  hazard            <- c()
  
  for (i in 1:n) {
    cumulative_hazard[i]  <- 0
    for (j in 1:J) {
      if (t[i] > interval_bounds[j] && t[i] <= interval_bounds[j+1]) {
        # Partial contribution for the current interval
        cumulative_hazard[i] <- cumulative_hazard[i] + lambda[j] * (t[i] - interval_bounds[j])
        current_interval     <- j
        break
      } else if (t[i] > interval_bounds[j+1]) {
        # Full contribution for earlier intervals
        cumulative_hazard[i] <- cumulative_hazard[i] + lambda[j] * (interval_bounds[j+1] - interval_bounds[j])
      }
    }
    hazard[i] <- lambda[current_interval]
  }
  return(list( h0 = hazard, H0 = cumulative_hazard ) )
}
#' @noRd
UpdateTau        <- function(survObj, priorPara, mcmcPara, ini){
  b.ini          <- ini$b.ini 
  sigmaSq        <- ini$sigmaSq
  etaSq          <- ini$etaSq	
  tauSq		       <- ini$tauSq
  nu.ind         <- NULL
  nu             <- sqrt(etaSq * sigmaSq *(b.ini^(-2)) )
  
  gam            <- c()
  for (j in 1:length(b.ini) ){
    repeat{
      gam[j]     <- rinvGauss(1, nu = nu[j], lambda = etaSq)
      if (gam[j] > 0) break
    }
    tauSq[j]  <- 1/gam[j]
  }
  
  return(tauSq)	
  
}
#' @noRd
UpdateTauStar    <- function( survObj, priorPara, mcmcPara, ini){
  beta.ini       <- ini$beta.ini
  sigmaStarSq    <- ini$sigmaStarSq
  etaStarSq      <- ini$etaStarSq	
  tauStarSq		   <- ini$tauStarSq
  nu.ind         <- NULL
  nu             <- sqrt(etaStarSq * sigmaStarSq *(beta.ini^(-2)) )
  gam            <- c()
  
  for (j in 1:length(beta.ini) ){
    repeat{
      gam[j]     <- rinvGauss(1, nu = nu[j], lambda = etaStarSq)
      if (gam[j] > 0) break
    }
    tauStarSq[j] <- 1/gam[j]
  }
  
  return(tauStarSq)	
  
}
#' @noRd
UpdateSigma      <- function( survObj, priorPara, mcmcPara, ini){
  b.ini          <- ini$b.ini 
  tauSq          <- ini$tauSq	
  sh.sig         <- length(b.ini)/2
  rate.sig       <- sum( (b.ini^2)/(2*tauSq) )
  sig.sq         <- rigamma(1, a = sh.sig, b = rate.sig)
  return(sig.sq)
}
#' @noRd
UpdateSigmaStar  <- function( survObj, priorPara, mcmcPara, ini){
  beta.ini       <- ini$beta.ini
  tauStarSq      <- ini$tauStarSq		
  sh.sig         <- length(beta.ini)/2
  rate.sig       <- sum( (beta.ini^2)/(2*tauStarSq) ) 
  sigStar.sq     <- rigamma(1, a = sh.sig, b = rate.sig)
  return(sigStar.sq)
}
#' @noRd
UpdateEta        <- function(survObj, priorPara, mcmcPara, ini){
  r1             <- priorPara$r1 
  delta1         <- priorPara$delta1 
  tauSq          <- ini$tauSq	
  etaSq	         <- rgamma(1, shape = r1 + length(tauSq), rate = delta1 + sum(tauSq)/2 )
  return(etaSq)
}
#' @noRd
UpdateEtaStar    <- function(survObj, priorPara, mcmcPara, ini){
  r2             <- priorPara$r2 
  delta2         <- priorPara$delta2 
  tauStarSq      <- ini$tauStarSq	
  etaStarSq	     <- rgamma(1, shape = r2 + length(tauStarSq), rate = delta2 + sum(tauStarSq)/2 )
  return(etaStarSq)
}