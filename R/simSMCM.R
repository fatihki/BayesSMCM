#' simSMFM
#' @title Function to generate data from Semi-parametric Mixture Cure Model (SMCM)
#'
#' @description
#' Simulate data from Semi-parametric Mixture Cure Model (SMCM) with the piecewise exponential baseline hazard function.
#'
#' @name simSMCM
#' 
#' @importFrom mvnfast rmvn 
#'
#' @param n Sample size
#' @param b.true True regression coefficients for the incidence part of the model.
#' @param beta.true True regression coefficients for the latency part of the model. 
#' @param baseline.hazard.rates Baseline hazard parameters for the piecewise exponential distribution.
#' @param intervals Vector of interval boundaries.
#' @param seed Random seed for reproducibility.
#' @param cens.start Start time for censoring distribution.
#' @param cens.end End time for censoring distribution.
#' @param prob.cov.X Vector of probabilities for generating binary covariates in \eqn{\bm X}, size at most length(beta.true).
#' @param prob.cov.Z Vector of probabilities for generating binary covariates in \eqn{\bm Z}, size at most length(b.true) - 1.
#' @param same.cov Logical indicating whether to use the same covariates for both \eqn{\bm X} and \eqn{\bm Z}.
#'
#' @return A list containing:
#' \item{b.true}{True regression coefficients for the incidence part of the model.}
#' \item{beta.true}{True regression coefficients for the latency part of the model.}
#' \item{baseline.hazard.rates}{Baseline hazard parameters for the piecewise exponential distribution.}
#' \item{X}{Covariate matrix for the incidence part of the model.}
#' \item{Z}{Covariate matrix for the latency part of the model.}
#' \item{observed_time}{Observed times to the event.}
#' \item{delta}{Censoring indicators for the event time (1=event occurred, 0=censored).}
#' \item{censor.rate}{Censoring rate for the simulated data.}
#' \item{intervals}{Vector of interval boundaries.}
#' 
#' @example
#' 
#' dat1 <- simSMCM(n = 300, b.true = c(0.4, 0.5, 0.1), beta.true = c(1, 0.2),
#'                 baseline.hazard.rates = c(0.1, 0.3), intervals = c(0.3, 0.6, 0.9), 
#'                 seed = 2026, cens.start = 30/365, cens.end = 32,
#'                 prob.cov.X = 0.5, prob.cov.Z = NULL, same.cov = TRUE)
#' 
#' @references 
#' Kızılaslan, F., Vitelli, V. (2025). Bayesian Semiparametric Mixture Cure (Frailty)
#' Model, arXiv \url{https://arxiv.org/abs/2512.08173}.
#'
#' @export
simSMCM <- function(n, b.true, beta.true, baseline.hazard.rates, intervals,
                       seed = NULL, cens.start, cens.end,
                       prob.cov.X, prob.cov.Z, same.cov = FALSE,
                       rho.X = 0, rho.Z = 0, cor.type = c("ar1", "cs", "block"), block.size = 10,
                       t.start = 0) {        
  
  cor.type <- match.arg(cor.type)
  if (!is.null(seed)) set.seed(seed)          
  
  ## ---- dimensions ----
  if (is.null(beta.true)) stop("Provide beta.true")
  p1 <- length(beta.true) 
  if (same.cov) {
    p2 <- p1 + 1
  } else {
    if (is.null(b.true)) stop("Provide b.true")
    p2 <- length(b.true)
  }
  if (same.cov && !is.null(b.true) && length(b.true) != p1 + 1)
    stop("same.cov = TRUE requires length(b.true) == length(beta.true) + 1.")
  
  cont.cov.X <- p1 - length(prob.cov.X)
  cont.cov.Z <- p2 - length(prob.cov.Z) - 1
  if (cont.cov.X < 0) stop("More binary X covariates than coefficients.")
  if (!same.cov && cont.cov.Z < 0) stop("More binary Z covariates than coefficients.")
  
  X  <- gen_cov(n, prob.cov.X, cont.cov.X, rho.X, cor.type, block.size)
  Zc <- if (same.cov) X else gen_cov(n, prob.cov.Z, cont.cov.Z, rho.Z, cor.type, block.size)
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  
  
  Z <- cbind(1, Zc)
  colnames(Z) <- c("Intercept", if (same.cov) colnames(X) else paste0("Z", seq_len(ncol(Zc))))
  
  pz.true <- drop(plogis(Z %*% b.true))         # inverse logit: P(susceptible)
  
  ## ---- survival times ----
  data    <- generate_piecewise_exponential_survival(cens.start, cens.end, intervals,
                                                     baseline.hazard.rates,
                                                     b_coef = b.true, betaa_coef = beta.true,
                                                     mX = X, mZ = Z, t.start = t.start )
  
  list(b.true = b.true, beta.true = beta.true,
       baseline.hazard.rates = baseline.hazard.rates,
       X = X, Z = Z,
       observed_time = data$observed_time, delta = data$delta,
       censor.rate = data$censor.rate,
       pz.true = pz.true, 
       cure.rate = 1 - mean(data$susceptible), susceptible = data$susceptible, 
       intervals = intervals)
}
#'
#' @noRd
generate_piecewise_exponential_survival <- function(cens.start, cens.end, intervals, baseline.hazard.rates, b_coef, betaa_coef, mX, mZ, t.start) {
  # n: Number of individuals
  # intervals: Vector of interval boundaries (e.g., c(0, 3, 6, 9))
  # baseline_hazard_rates: Baseline hazard rates for each interval (e.g., c(0.2, 0.1, 0.05))
  # beta: Vector of regression coefficients (matching the number of covariates)
  # covariates: Matrix of covariates (n x p), where p is the number of covariates
  n                 <- nrow(mX)
  t_star            <- rep(NA_real_, n)
  pi_Z              <- rep(NA_real_, n)
  susceptible <- rep(NA_integer_, n)
  
  if (length(intervals) != (length(baseline.hazard.rates) + 1)) {
    stop("The length of intervals should be one more than the length of baseline.hazard.rates.")
  }
  
  if (ncol(mX) != length(betaa_coef)) {
    stop("The length of beta should match the number of covariate columns.")
  }
  
  survival_times     <- rep(NA_real_, n)
  covariate_effects  <- mX %*% betaa_coef 
  
  for (i in 1:n) {
    # Calculate individual-specific hazard rates based on covariates
    lambda_i         <- baseline.hazard.rates * exp(covariate_effects[i])  # Adjusted hazard rates for individual i
    T                <- t.start  # Initialize survival time for individual i
    
    # Simulate survival time across the piecewise intervals
    for (j in 1:length(lambda_i)) {
      # Calculate the survival time in the current interval
      u              <- runif(1)
      t_j            <- -log(u) / lambda_i[j]  # Exponential survival time in the current interval
      
      # Check if the event occurs within the current interval
      if (j == length(lambda_i)  || T + t_j    < intervals[j + 1]) {
        T            <- T + t_j
        break
      } else {
        # If event does not occur in this interval, move to the next interval
        T             <- intervals[j + 1]
      }
    }
    survival_times[i] <- T
  }
  
  for (i in 1:n) {
    pi_Z[i]        <- plogis(sum(mZ[i, ] * b_coef))           # P(susceptible)
    susceptible[i] <- as.integer(runif(1) < pi_Z[i])
    t_star[i]      <- if (susceptible[i] == 1L) survival_times[i] else Inf
  }

  censor_time     <-runif(n, cens.start, cens.end)
  observed_time   <- pmin(t_star, censor_time)
  event          <- as.numeric(t_star <= censor_time)  # 1 = event observed, 0 = censored
  
  list(observed_time = observed_time, delta = event, censor.rate = 1 - (sum(event) / n),
       susceptible = susceptible)
  
}
#'
#' Correlation matrix for the continuous covariates
#'   "ar1"   : rho^|i-j|
#'   "cs"    : compound symmetry (all pairs rho)
#'   "block" : rho within blocks of size block.size, 0 between blocks
#' @examples
#'make_cor(p=6, rho=0.5, type = "ar1" )
# make_cor(p=6, rho=0.5, type = "cs" )
# make_cor(p=6, rho=0.5, type = "block", block.size = 2 )
#' @noRd
make_cor <- function(p, rho, type = "ar1", block.size = 10) {
  if (rho == 0) return(diag(p))
  switch(type,
         ar1   = rho^abs(outer(seq_len(p), seq_len(p), "-")),
         cs    = { S <- matrix(rho, p, p); diag(S) <- 1; S },
         block = { blk <- (seq_len(p) - 1) %/% block.size
         S <- outer(blk, blk, "==") * rho; diag(S) <- 1; S })
}
#'
#'
## ---- covariates: binary columns first, then (correlated) N(0, Sigma) columns ----
#' @example 
#' gen_cov(n=10, probs=0.5, n.cont=2, rho = 0, cor.type = "ar1", block.size = 10)
#' @noRd
gen_cov <- function(n, probs, n.cont, rho = 0, cor.type = "ar1", block.size = 10) {
  B <- if (length(probs) > 0) matrix(sapply(probs, function(x) rbinom(n, 1, x)), nrow = n) else NULL
  C <- if (n.cont > 0) mvnfast::rmvn(n, mu = rep(0, n.cont),
                            sigma = make_cor(n.cont, rho, cor.type, block.size)) else NULL
  cbind(B, C)
}
