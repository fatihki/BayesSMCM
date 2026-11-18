#' simSMCFM
#' @title Data generation from the mixture cure FRAILTY model (SMCFM)
#'
#' @description
#' Same interface as \code{simSMCM}, plus the frailty parameter \code{theta.true}.
#' Susceptible (uncured) subjects, with probability \eqn{\pi(z)=\mathrm{plogis}(z^\top b)}, have
#' conditional hazard \eqn{w_i h_0(t) e^{x_i^\top\beta}} with a piecewise-exponential \eqn{h_0} and
#' frailty \eqn{w_i \sim \mathrm{Gamma}(\theta, \theta)} (mean 1, variance \eqn{1/\theta}), so that
#' marginally \eqn{S(t|x) = (1 + H_0(t)e^{x^\top\beta}/\theta)^{-\theta}}.
#' Cured subjects never experience the event. Censoring is Uniform(cens.start, cens.end).
#'
#' @param theta.true Frailty parameter theta > 0 (variance of the frailty = 1 / theta).
#' @param intervals Interval boundaries of the baseline hazard (length J + 1);
#'   \code{baseline.hazard.rates} has length J. The last rate is extended beyond the last boundary.
#' @inheritParams simSMCM
#'
#' @return A list with b.true, beta.true, theta.true, baseline.hazard.rates, X, Z, observed_time,
#'   delta, censor.rate, pz.true, cure.rate (realised), intervals, and the latent
#'   susceptible indicator and frailty values (\code{susceptible}, \code{frailty}).
#'   
#' @example 
#' dat1 <- simSMCFM(n = 300, b.true = c(0.4, 0.5, 0.1), beta.true = c(1, 0.2), theta.true = 1,
#'                  baseline.hazard.rates = c(0.1, 0.3), intervals = c(0.3, 0.6, 0.9), 
#'                  seed = 2026, cens.start = 30/365, cens.end = 32,
#'                  prob.cov.X = 0.5, prob.cov.Z = NULL, same.cov = TRUE)
#' str(dat1)
#' 
#' @references 
#' Kızılaslan, F., Vitelli, V. (2025). Bayesian Semiparametric Mixture Cure (Frailty)
#' Model, arXiv \url{https://arxiv.org/abs/2512.08173}.
#' 
#' @export
simSMCFM <- function(n, b.true, beta.true, theta.true, baseline.hazard.rates, intervals,
                        seed = NULL, cens.start, cens.end,
                        prob.cov.X, prob.cov.Z, same.cov = FALSE,
                        rho.X = 0, rho.Z = 0, cor.type = c("ar1", "cs", "block"), block.size = 10,
                        t.start = 0) {
  
  cor.type <- match.arg(cor.type)
  if (!is.null(seed)) set.seed(seed)
  
  ## ---- checks ----
  if (!is.numeric(theta.true) || length(theta.true) != 1L || theta.true <= 0)
    stop("theta.true must be a single positive number.")
  if (length(intervals) != length(baseline.hazard.rates) + 1L)
    stop("length(intervals) must equal length(baseline.hazard.rates) + 1.")
  if (any(baseline.hazard.rates <= 0) || any(diff(intervals) <= 0))
    stop("baseline.hazard.rates must be positive and intervals strictly increasing.")
  
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
  
  pz.true <- drop(plogis(Z %*% b.true))            # P(susceptible)
  
  ## ---- survival times ----
  data <- generate_piecewise_exponential_frailty_survival(
    cens.start, cens.end, intervals, baseline.hazard.rates,
    b_coef = b.true, betaa_coef = beta.true, theta = theta.true,
    mX = X, mZ = Z, t.start = t.start)
  
  list(b.true = b.true, beta.true = beta.true, theta.true = theta.true,
       baseline.hazard.rates = baseline.hazard.rates,
       X = X, Z = Z,
       observed_time = data$observed_time, delta = data$delta,
       censor.rate = data$censor.rate,
       pz.true = pz.true, cure.rate = 1 - mean(data$susceptible),
       susceptible = data$susceptible, frailty = data$frailty,
       intervals = intervals)
}
#'
#' Piecewise-exponential survival times with Gamma frailty and a cured fraction
#'
#' Same construction as \code{generate_piecewise_exponential_survival} (SMCM), with one change:
#' every subject gets a frailty w_i ~ Gamma(theta, theta) (mean 1, variance 1/theta) that multiplies
#' the whole hazard, i.e. conditional on w_i the hazard is w_i * h0(t) * exp(x_i'beta).
#' Marginally (frailty integrated out) S(t|x) = (1 + H0(t) exp(x'beta) / theta)^(-theta).
#'
#' @param theta Frailty parameter (> 0).
#' @param t.start Time at which the piecewise simulation starts (default 0, i.e. event times follow the
#'   model exactly). A positive value acts as a minimum follow-up offset (the SMCM generator used 30/365).
#' @noRd
generate_piecewise_exponential_frailty_survival <- function(cens.start, cens.end, intervals,
                                                            baseline.hazard.rates,
                                                            b_coef, betaa_coef, theta, mX, mZ,
                                                            t.start = 0) {
  n          <- nrow(mX)
  t_star     <- rep(NA_real_, n)
  pi_Z       <- rep(NA_real_, n)
  frailty    <- rep(NA_real_, n)
  susceptible <- rep(NA_integer_, n)
  
  if (length(intervals) != (length(baseline.hazard.rates) + 1)) {
    stop("The length of intervals should be one more than the length of baseline.hazard.rates.")
  }
  if (ncol(mX) != length(betaa_coef)) {
    stop("The length of beta should match the number of covariate columns.")
  }
  if (!is.numeric(theta) || length(theta) != 1L || theta <= 0) stop("theta must be a positive number.")
  
  survival_times    <- rep(NA_real_, n)
  covariate_effects <- mX %*% betaa_coef                     # n x 1 vector of covariate effects
  J                 <- length(baseline.hazard.rates)
  
  for (i in 1:n) {
    # frailty: multiplies the whole hazard of subject i
    frailty[i] <- rgamma(1, shape = theta, rate = theta)
    lambda_i   <- baseline.hazard.rates * exp(covariate_effects[i]) * frailty[i]
    Ti         <- t.start
    
    # piecewise exponential simulation (memoryless property within each interval)
    for (j in 1:J) {
      if (j < J && intervals[j + 1] <= Ti) next               # start time already beyond this interval
      t_j <- -log(runif(1)) / lambda_i[j]
      if (j == J || Ti + t_j < intervals[j + 1]) {
        Ti <- Ti + t_j
        break
      } else {
        Ti <- intervals[j + 1]                                # no event here: move to next interval
      }
    }
    survival_times[i] <- Ti
  }
  
  for (i in 1:n) {
    pi_Z[i]        <- plogis(sum(mZ[i, ] * b_coef))           # P(susceptible)
    susceptible[i] <- as.integer(runif(1) < pi_Z[i])
    t_star[i]      <- if (susceptible[i] == 1L) survival_times[i] else Inf
  }
  
  censor_time   <- runif(n, cens.start, cens.end)
  observed_time <- pmin(t_star, censor_time)
  event         <- as.numeric(t_star <= censor_time)                 # 1 = event observed, 0 = censored
  
  list(observed_time = observed_time, delta = event, censor.rate = 1 - (sum(event) / n),
       susceptible = susceptible, frailty = frailty)
}
