#' simWMCM
#' @title Data generation from the parametric Weibull mixture cure model (WMCM)
#'
#' @description
#' Same interface as \code{simSMCM.R1} / \code{simSMCFM.R1}, but the latency part is a parametric
#' Weibull proportional-hazards model. Susceptible (uncured) subjects, with probability
#' \eqn{\pi(z)=\mathrm{plogis}(z^\top b)}, have survival function
#' \deqn{S(t|x) = \exp\{-\mathrm{scale}\; t^{\mathrm{shape}}\; e^{x^\top\beta}\},}
#' i.e. baseline hazard \eqn{h_0(t) = \mathrm{scale}\cdot\mathrm{shape}\cdot t^{\mathrm{shape}-1}} and baseline
#' cumulative hazard \eqn{H_0(t) = \mathrm{scale}\cdot t^{\mathrm{shape}}}. Cured subjects never experience the
#' event. Censoring is Uniform(cens.start, cens.end), optionally combined with an exponential loss to follow-up.
#'
#' Covariates are generated exactly as in \code{simSMCM.R1} (binary columns first, then continuous
#' N(0, Sigma) columns with AR(1), compound-symmetry or block correlation), using \code{gen_cov()} and
#' \code{make_cor()} from simSMCM_R1.R (they are NOT redefined here), and the package 'mvnfast'.
#'
#' @param n Sample size.
#' @param b.true,beta.true True coefficients of the incidence (length p2, first = intercept) and latency
#'   (length p1) parts.
#' @param surv_times,surv_val The Weibull baseline is defined by its survival curve at one or two time points,
#'   \eqn{S_0(t_k)} = \code{surv_val[k]} at \code{surv_times[k]}. One point gives an exponential baseline
#'   (shape = 1); two points (t1 < t2, S0(t1) > S0(t2)) give a Weibull baseline.
#' @param weibull.scale,weibull.shape Alternative: give the Weibull scale (rate-type, H0(t) = scale * t^shape)
#'   and shape directly. They take precedence over \code{surv_times}/\code{surv_val}.
#' @param seed Random seed (set once at the start, as in the other generators).
#' @param cens.start,cens.end Censoring times are Uniform(cens.start, cens.end).
#' @param prob.cov.X,prob.cov.Z Probabilities of the binary covariates in X and Z (the remaining columns
#'   are continuous).
#' @param same.cov If TRUE, Z = (1, X); requires \code{length(b.true) == length(beta.true) + 1}.
#' @param rho.X,rho.Z,cor.type,block.size Correlation among the continuous covariates (see \code{make_cor}).
#' @param t.start Time at which the hazard process starts (default 0: event times follow the model exactly).
#'   A positive value acts as a minimum follow-up offset: event times are drawn conditionally on survival
#'   beyond \code{t.start} (the same idea as in the piecewise generators). It replaces the old ad hoc rule
#'   "add 30/365 to times shorter than 30/365".
#' @param loss.surv Optional number in (0, 1). If given, an additional exponential loss to follow-up is
#'   applied with P(no loss by cens.end) = loss.surv (the old simWMCM used 0.8). Default NULL: uniform
#'   censoring only, as in simSMCM.R1 and simSMCFM.R1.
#'
#' @return A list with b.true, beta.true, weibull (list with scale and shape), surv_times, surv_val, X, Z,
#'   observed_time, delta, censor.rate (realised), pz.true, cure.rate (realised), susceptible (latent
#'   indicator), t.start, loss.surv.
#'
#' @examples
#' # Low dimension, same covariates in both parts, increasing Weibull hazard (S0(5) = 0.9, S0(15) = 0.5)
#' dat1 <- simWMCM(n = 300, b.true = c(0.4, 0.5, 0.1), beta.true = c(1, 0.2),
#'                    surv_times = c(5, 15), surv_val = c(0.9, 0.5), seed = 2026,
#'                    cens.start = 30/365, cens.end = 32,
#'                    prob.cov.X = 0.5, prob.cov.Z = NULL, same.cov = TRUE)
#'
#' # High dimension, sparse truth, AR(1)-correlated continuous covariates (rho = 0.3)
#' p <- 50
#' dat2 <- simWMCM(n = 500, b.true = c(0.4, 1, -1, rep(0, p - 2)),
#'                    beta.true = c(1, -0.5, 0.5, rep(0, p - 3)),
#'                    weibull.scale = 0.01, weibull.shape = 1.5, seed = 2026,
#'                    cens.start = 30/365, cens.end = 32,
#'                    prob.cov.X = rep(0.5, 3), prob.cov.Z = NULL, same.cov = TRUE,
#'                    rho.X = 0.3, cor.type = "ar1")
#'                    
#' @references 
#' Kızılaslan, F., Vitelli, V. (2025). Bayesian Semiparametric Mixture Cure (Frailty)
#' Model, arXiv \url{https://arxiv.org/abs/2512.08173}.
#'                    
#' @export
simWMCM <- function(n, b.true, beta.true,
                       surv_times = NULL, surv_val = NULL,
                       weibull.scale = NULL, weibull.shape = NULL,
                       seed = NULL, cens.start, cens.end,
                       prob.cov.X, prob.cov.Z, same.cov = FALSE,
                       rho.X = 0, rho.Z = 0, cor.type = c("ar1", "cs", "block"), block.size = 10,
                       t.start = 0, loss.surv = NULL) {
  
  cor.type <- match.arg(cor.type)
  if (!is.null(seed)) set.seed(seed)
  
  ## ---- checks ----
  if (is.null(beta.true)) stop("Provide beta.true")
  if (is.null(b.true))    stop("Provide b.true")
  if (!is.numeric(t.start) || length(t.start) != 1L || t.start < 0)
    stop("t.start must be a single non-negative number.")
  if (cens.start < 0 || cens.end <= cens.start)
    stop("Need 0 <= cens.start < cens.end.")
  if (!is.null(loss.surv) &&
      (!is.numeric(loss.surv) || length(loss.surv) != 1L || loss.surv <= 0 || loss.surv >= 1))
    stop("loss.surv must be NULL or a single number in (0, 1).")
  
  wb <- weibull_params(surv_times, surv_val, weibull.scale, weibull.shape)
  
  ## ---- dimensions (identical to simSMCM.R1) ----
  p1 <- length(beta.true)
  p2 <- if (same.cov) p1 + 1 else length(b.true)
  if (same.cov && length(b.true) != p1 + 1)
    stop("same.cov = TRUE requires length(b.true) == length(beta.true) + 1.")
  
  cont.cov.X <- p1 - length(prob.cov.X)
  cont.cov.Z <- p2 - length(prob.cov.Z) - 1
  if (cont.cov.X < 0) stop("More binary X covariates than coefficients.")
  if (!same.cov && cont.cov.Z < 0) stop("More binary Z covariates than coefficients.")
  
  ## ---- covariates (binary first, then correlated N(0, Sigma) columns) ----
  X  <- gen_cov(n, prob.cov.X, cont.cov.X, rho.X, cor.type, block.size)
  Zc <- if (same.cov) X else gen_cov(n, prob.cov.Z, cont.cov.Z, rho.Z, cor.type, block.size)
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  
  Z <- cbind(1, Zc)
  colnames(Z) <- c("Intercept", if (same.cov) colnames(X) else paste0("Z", seq_len(ncol(Zc))))
  
  pz.true <- drop(plogis(Z %*% b.true))                    # P(susceptible)
  
  ## ---- survival data ----
  data <- generate_weibull_mixture_cure_survival(
    cens.start, cens.end, wb$scale, wb$shape,
    b_coef = b.true, betaa_coef = beta.true,
    mX = X, mZ = Z, t.start = t.start, loss.surv = loss.surv)
  
  list(b.true = b.true, beta.true = beta.true,
       weibull = wb, surv_times = surv_times, surv_val = surv_val,
       X = X, Z = Z,
       observed_time = data$observed_time, delta = data$delta,
       censor.rate = data$censor.rate,
       pz.true = pz.true, cure.rate = 1 - mean(data$susceptible),
       susceptible = data$susceptible,
       t.start = t.start, loss.surv = loss.surv)
}


#' Weibull parameters from two (or one) points of the baseline survival curve, or given directly
#'
#' H0(t) = scale * t^shape. From S0(t1) = s1 and S0(t2) = s2:
#'   shape = log(z2 / z1) / log(t2 / t1),  scale = z1 / t1^shape,  with z_k = -log(s_k).
#' One point: exponential (shape = 1, scale = z1 / t1).
#' @noRd
weibull_params <- function(surv_times, surv_val, scale = NULL, shape = NULL) {
  
  if (!is.null(scale) || !is.null(shape)) {
    if (is.null(scale) || is.null(shape))
      stop("Provide both weibull.scale and weibull.shape.")
    if (!is.numeric(scale) || !is.numeric(shape) || length(scale) != 1L || length(shape) != 1L ||
        scale <= 0 || shape <= 0)
      stop("weibull.scale and weibull.shape must be single positive numbers.")
    return(list(scale = scale, shape = shape))
  }
  
  if (is.null(surv_times) || is.null(surv_val))
    stop("Provide either (surv_times, surv_val) or (weibull.scale, weibull.shape).")
  if (length(surv_times) != length(surv_val) || !(length(surv_times) %in% 1:2))
    stop("surv_times and surv_val must have the same length, 1 or 2.")
  if (any(surv_times <= 0) || any(surv_val <= 0 | surv_val >= 1))
    stop("surv_times must be positive and surv_val must lie in (0, 1).")
  
  z <- -log(surv_val)
  if (length(z) == 1L) {
    shape <- 1
  } else {
    if (surv_times[2] <= surv_times[1] || surv_val[2] >= surv_val[1])
      stop("Need surv_times[1] < surv_times[2] and surv_val[1] > surv_val[2].")
    shape <- log(z[2] / z[1]) / log(surv_times[2] / surv_times[1])
  }
  list(scale = z[1] / surv_times[1]^shape, shape = shape)
}
#'
#' Weibull event times, cure status and censoring
#'
#' Inverse-transform sampling from S(t|x) = exp(-scale t^shape e^{x'beta}), conditional on survival
#' beyond t.start (t.start = 0: unconditional). With E = -log(U) ~ Exp(1):
#'   scale * (t^shape - t.start^shape) * e^{x'beta} = E
#'   => t = (t.start^shape + E * exp(-x'beta) / scale)^(1 / shape).
#' @noRd
generate_weibull_mixture_cure_survival <- function(cens.start, cens.end, wb.scale, wb.shape,
                                                   b_coef, betaa_coef, mX, mZ,
                                                   t.start = 0, loss.surv = NULL) {
  n <- nrow(mX)
  if (ncol(mX) != length(betaa_coef))
    stop("The length of beta should match the number of covariate columns.")
  if (ncol(mZ) != length(b_coef))
    stop("The length of b should match the number of columns of Z.")

  eta     <- drop(mX %*% betaa_coef)                       # latency linear predictor
  E       <- -log(runif(n))
  t.event <- (t.start^wb.shape + E * exp(-eta) / wb.scale)^(1 / wb.shape)
  if (anyNA(t.event) || any(t.event < 0))
    stop("Invalid event times generated; check the Weibull parameters and the coefficients.")
  
  pi_Z        <- plogis(drop(mZ %*% b_coef))               # P(susceptible)
  susceptible <- as.integer(runif(n) < pi_Z)
  t_star      <- ifelse(susceptible == 1L, t.event, Inf)   # cured: never fail
  
  censor_time <- runif(n, cens.start, cens.end)
  # if (!is.null(loss.surv)) {                               # optional exponential loss to follow-up
  #   loss.rate   <- -log(loss.surv) / cens.end              # P(no loss by cens.end) = loss.surv
  #   censor_time <- pmin(censor_time, rexp(n, rate = loss.rate))
  # }
  observed_time <- pmin(t_star, censor_time)
  event         <- as.numeric(t_star <= censor_time)       # 1 = event observed, 0 = censored
  
  list(observed_time = observed_time, delta = event,
       censor.rate = 1 - (sum(event) / n), susceptible = susceptible)
}
