#' @title Logarithm of the likelihood for SMCFM
#'
#' @description
#' Compute the logarithm of the likelihood of SMCFM for the right-censored observed survival data. 
#'
#' @name logLik_SMCFM
#'
#' @param survObj An object of class SMCFM containing survival data and covariates.
#' @param priorPara A list of prior parameters including interval bounds.
#' @param b.coef regression coefficients for the incidence part of the model.
#' @param beta.coef regression coefficients for the latency part of the model.
#' @param lambda.coef baseline hazard parameters for the piecewise exponential distribution.
#' @param theta.coef gamma distribution parameter for the frailty distribution.
#' 
#' @return Numeric value of the log-likelihood.
#' 
#' @export
logLik_SMCFM <- function(survObj, priorPara, b.coef, beta.coef, lambda.coef, theta.coef ){
  Z          <- survObj$Z
  X          <- survObj$X
  y          <- survObj$t
  delta      <- survObj$di
  s          <- priorPara$s
  # Compute cumulative hazard and hazard for each t[i]
  H0_h0      <- cumulative_hazard(y, interval_bounds = s, lambda = lambda.coef)
  H0         <- H0_h0$H0 
  h0         <- H0_h0$h0 #
  Zb         <- Z %*% b.coef 
  exp.Zb     <- exp(Zb)
  Xbeta      <- X %*% beta.coef
  exp.Xbeta  <- exp(Xbeta)
  H1         <- H0 * exp.Xbeta  
  H20        <- H1 / theta.coef  
  H2         <- 1 + H20          
  term1      <- ( Zb + log(h0) + Xbeta - (theta.coef+1) * log1p(H20) ) * delta 
  term2      <- log1p(exp.Zb * H2^(-theta.coef)) * (1-delta) 
  logL       <- term1 + term2 - log1p(exp.Zb)                                    # log-likelihood for each observation
  logL       <- t(logL)
  colnames(logL) <- paste0("V", 1:length(logL))
  
  return(logL)
  
}