# Optimized implementation for the Bayesian semiparametric
# mixture cure FRAILTY model (SMCFM).
#
# Model: S_pop(t | x, z) = 1 - pi(z) + pi(z) * L( H0(t) exp(x'beta) ),
#   pi(z) = plogis(z'b)                                   (incidence, coefficients b)
#   L(s)  = (1 + s / theta)^(-theta)                      (Gamma frailty, mean 1, var 1/theta,
#                                                          acting on the latency part)
#   H0    = piecewise-exponential cumulative baseline hazard.
#
# log S_frailty(t_i) = -theta * log(1 + H1_i / theta)
smcfm_logS <- function(H1, theta) -theta * log1p(H1 / theta)

# b-dependent part. `logS` = smcfm_logS(H1, theta) does not depend on b.
smcfm_loglik_b <- function(Zb, logS, event) {
  cens <- !event
  out <- -sum(softplus(Zb)) + sum(Zb[event])
  if (any(cens)) out <- out + sum(softplus(Zb[cens] + logS[cens]))
  out
}

# beta-dependent part.
smcfm_loglik_beta <- function(Xbeta, Zb, H0, theta, event) {
  l1p  <- log1p(H0 * exp(Xbeta) / theta)
  cens <- !event
  out  <- sum(Xbeta[event] - (theta + 1) * l1p[event]) +
    sum(softplus(Zb[cens] - theta * l1p[cens]))
  if (!is.finite(out)) -Inf else out
}

# lambda-dependent part.
smcfm_loglik_lambda <- function(H0, h0, Xbeta, Zb, theta, event) {
  if (any(h0 <= 0) || any(!is.finite(h0))) return(-Inf)
  l1p  <- log1p(H0 * exp(Xbeta) / theta)
  cens <- !event
  out  <- sum(log(h0[event]) - (theta + 1) * l1p[event]) +
    sum(softplus(Zb[cens] - theta * l1p[cens]))
  if (!is.finite(out)) -Inf else out
}

# theta-dependent part (H1 = H0 * exp(Xbeta) is fixed).
smcfm_loglik_theta <- function(theta, H1, Zb, event) {
  l1p  <- log1p(H1 / theta)
  cens <- !event
  out  <- -(theta + 1) * sum(l1p[event]) + sum(softplus(Zb[cens] - theta * l1p[cens]))
  if (!is.finite(out)) -Inf else out
}

# Full per-observation log-likelihood (1 x n matrix), for LOO / DIC / LPML and validation.
logLik_SMCFM_fast <- function(survObj, priorPara, b.coef, beta.coef, lambda.coef, theta.coef,
                              hazard_info = NULL) {
  if (is.null(hazard_info)) {
    hazard_info <- cumulative_hazard_matrix(t = survObj$t, interval_bounds = priorPara$s,
                                            lambda = lambda.coef)
  }
  Zb    <- as.vector(survObj$Z %*% b.coef)
  Xbeta <- as.vector(survObj$X %*% beta.coef)
  hz    <- smcm_hazard_from_state(lambda.coef, hazard_info)
  l1p   <- log1p(hz$H0 * exp(Xbeta) / theta.coef)

  event <- survObj$di == 1
  cens  <- !event
  ll <- -softplus(Zb)
  if (any(event)) {
    ll[event] <- ll[event] + Zb[event] + log(hz$h0[event]) + Xbeta[event] -
      (theta.coef + 1) * l1p[event]
  }
  if (any(cens)) {
    ll[cens] <- ll[cens] + softplus(Zb[cens] - theta.coef * l1p[cens])
  }
  matrix(ll, nrow = 1L, dimnames = list(NULL, paste0("V", seq_along(ll))))
}

# Backward-compatible name used by mcmc.dic.lpml() and older code.
logLik_SMCFM <- function(survObj, priorPara, b.coef, beta.coef, lambda.coef, theta.coef) {
  logLik_SMCFM_fast(survObj, priorPara, b.coef, beta.coef, lambda.coef, theta.coef)
}

##=============================================================================##
##===================== Metropolis-Hastings updates============================##
##=============================================================================##

UpdateRP.b.SMCFM.fast <- function(Z, b, Zb, logS, event, sd.b,
                                  proposal.sd = rep(1, length(b))) {
  p <- length(b)
  if (length(proposal.sd) == 1L) proposal.sd <- rep(proposal.sd, p)
  accept <- integer(p)
  current.loglik <- smcfm_loglik_b(Zb, logS, event)

  for (k in seq_len(p)) {
    b.prop  <- rnorm(1L, mean = b[k], sd = proposal.sd[k])
    Zb.prop <- Zb + Z[, k] * (b.prop - b[k])
    prop.loglik <- smcfm_loglik_b(Zb.prop, logS, event)

    logR <- prop.loglik - current.loglik +
      dnorm(b.prop, 0, sd.b[k], log = TRUE) - dnorm(b[k], 0, sd.b[k], log = TRUE)

    if (!is.na(logR) && log(runif(1L)) < logR) {
      b[k] <- b.prop;  Zb <- Zb.prop
      current.loglik <- prop.loglik;  accept[k] <- 1L
    }
  }
  list(b = b, Zb = Zb, loglik = current.loglik, accept = accept)
}

UpdateRP.beta.SMCFM.fast <- function(X, beta, Xbeta, Zb, H0, theta, event, sd.beta,
                                     proposal.sd = rep(1, length(beta))) {
  p <- length(beta)
  if (length(proposal.sd) == 1L) proposal.sd <- rep(proposal.sd, p)
  accept <- integer(p)
  current.loglik <- smcfm_loglik_beta(Xbeta, Zb, H0, theta, event)

  for (k in seq_len(p)) {
    beta.prop  <- rnorm(1L, mean = beta[k], sd = proposal.sd[k])
    Xbeta.prop <- Xbeta + X[, k] * (beta.prop - beta[k])
    prop.loglik <- smcfm_loglik_beta(Xbeta.prop, Zb, H0, theta, event)

    logR <- prop.loglik - current.loglik +
      dnorm(beta.prop, 0, sd.beta[k], log = TRUE) - dnorm(beta[k], 0, sd.beta[k], log = TRUE)

    if (!is.na(logR) && log(runif(1L)) < logR) {
      beta[k] <- beta.prop;  Xbeta <- Xbeta.prop
      current.loglik <- prop.loglik;  accept[k] <- 1L
    }
  }
  list(beta = beta, Xbeta = Xbeta, loglik = current.loglik, accept = accept)
}

# lambda_k: Gaussian random walk on log(lambda_k); Jacobian log(lambda*) - log(lambda) in the ratio.
UpdateBH.lambda.SMCFM.fast <- function(lambda, H0, h0, hazard_info, Xbeta, Zb, theta, event,
                                       a, b, proposal.sd = rep(1, length(lambda))) {
  J <- length(lambda)
  if (length(proposal.sd) == 1L) proposal.sd <- rep(proposal.sd, J)
  if (length(proposal.sd) != J) stop("proposal.sd must have length 1 or length(lambda).")
  if (any(!is.finite(proposal.sd)) || any(proposal.sd <= 0)) {
    stop("proposal.sd must contain positive finite values.")
  }
  if (any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("lambda must contain positive finite values.")
  }

  accept <- integer(J)
  current.loglik <- smcfm_loglik_lambda(H0, h0, Xbeta, Zb, theta, event)

  for (k in seq_len(J)) {
    lambda.current <- lambda[k]
    log.lambda.current <- log(lambda.current)
    log.lambda.prop <- log.lambda.current + rnorm(1L, 0, proposal.sd[k])
    if (!is.finite(log.lambda.prop)) next
    lambda.prop <- exp(log.lambda.prop)
    if (!is.finite(lambda.prop) || lambda.prop <= 0) next

    # incremental update: H0* = H0 + A[, k] (lambda_k* - lambda_k)
    H0.prop <- H0 + hazard_info$A[, k] * (lambda.prop - lambda.current)
    if (any(!is.finite(H0.prop)) || any(H0.prop < 0)) next

    h0.prop <- h0
    idx <- hazard_info$interval_id == k
    if (any(idx)) h0.prop[idx] <- lambda.prop

    prop.loglik <- smcfm_loglik_lambda(H0.prop, h0.prop, Xbeta, Zb, theta, event)
    if (!is.finite(prop.loglik)) next

    logR <- prop.loglik - current.loglik +
      dgamma(lambda.prop, shape = a, rate = b, log = TRUE) -
      dgamma(lambda.current, shape = a, rate = b, log = TRUE) +
      log.lambda.prop - log.lambda.current

    if (!is.na(logR) && log(runif(1L)) < logR) {
      lambda[k] <- lambda.prop;  H0 <- H0.prop;  h0 <- h0.prop
      current.loglik <- prop.loglik;  accept[k] <- 1L
    }
  }
  list(lambda = lambda, H0 = H0, h0 = h0, loglik = current.loglik, accept = accept)
}

# theta: Gaussian random walk on log(theta), prior theta ~ Gamma(shape = c, rate = d).
UpdateTheta.SMCFM.fast <- function(theta, H1, Zb, event, c.shape, d.rate, proposal.sd = 0.5) {
  current.loglik <- smcfm_loglik_theta(theta, H1, Zb, event)

  log.theta.prop <- log(theta) + rnorm(1L, 0, proposal.sd)
  theta.prop     <- exp(log.theta.prop)
  if (!is.finite(theta.prop) || theta.prop <= 0) {
    return(list(theta = theta, loglik = current.loglik, accept = 0L))
  }
  prop.loglik <- smcfm_loglik_theta(theta.prop, H1, Zb, event)

  logR <- prop.loglik - current.loglik +
    dgamma(theta.prop, shape = c.shape, rate = d.rate, log = TRUE) -
    dgamma(theta,      shape = c.shape, rate = d.rate, log = TRUE) +
    log.theta.prop - log(theta)

  if (!is.na(logR) && log(runif(1L)) < logR) {
    return(list(theta = theta.prop, loglik = prop.loglik, accept = 1L))
  }
  list(theta = theta, loglik = current.loglik, accept = 0L)
}
#'
#'
#'
#' @title Function to MCMC algorithm for Semi-parametric Mixture Cure Frailty Model (SMCFM)
#'
#' @description
#' This function obtains posterior samples for all model parameters using an MCMC algorithm 
#' for a single chain when applying the \code{MCMC.SMCFM.fast} framework  in
#'  Kızılaslan and Vitelli (2025) to the data. 
#'
#' @name MCMC.SMCFM.fast
#' @export
MCMC.SMCFM.fast <- function(survObj, priorPara, initial,
                            nIter = 6000L, warmup = 1000L, thin = 1L,
                            save_loglik = 0L, store_tau = TRUE,
                            proposal_sd_b = 1,
                            proposal_sd_beta = 1,
                            proposal_sd_lambda = 0.5,
                            proposal_sd_theta = 0.5,
                            adapt = TRUE, target_accept = 0.35,
                            refresh_every = 200L,
                            progress = interactive()) {
  stopifnot(nIter >= 1L, warmup >= 0L, warmup < nIter, thin >= 1L)

  X <- as.matrix(survObj$X)
  Z <- as.matrix(survObj$Z)
  delta <- as.integer(survObj$di)
  event <- delta == 1L
  n <- nrow(X)
  p_beta <- ncol(X)
  p_b <- ncol(Z)
  J <- priorPara$J
  if (length(priorPara$s) - 1L != J) stop("priorPara$J does not match length(priorPara$s)-1.")
  if (is.null(priorPara$c) || is.null(priorPara$d)) {
    stop("priorPara must contain c and d (Gamma prior of the frailty parameter theta).")
  }

  beta <- as.numeric(initial$beta.ini)
  bcoef <- as.numeric(initial$b.ini)
  lambda <- as.numeric(initial$lambda.ini)
  theta <- as.numeric(initial$theta.ini)
  sigmaSq <- initial$sigmaSq
  sigmaStarSq <- initial$sigmaStarSq
  etaSq <- initial$etaSq
  etaStarSq <- initial$etaStarSq
  tauSq <- as.numeric(initial$tauSq)
  tauStarSq <- as.numeric(initial$tauStarSq)

  ## N(0, intercept_sd^2) prior and is excluded from the lasso hierarchy unless
  ## priorPara$penalize_intercept is TRUE.
  has_int      <- p_b > 0L && all(Z[, 1L] == 1)
  unpen_b      <- if (!isTRUE(priorPara$penalize_intercept) && has_int) 1L else integer(0)
  pen_b        <- setdiff(seq_len(p_b), unpen_b)
  intercept_sd <- if (is.null(priorPara$intercept_sd)) 10 else priorPara$intercept_sd
  tauSq[unpen_b] <- 1
  a_sigma <- if (is.null(priorPara$a_sigma)) 3 else priorPara$a_sigma
  b_sigma <- if (is.null(priorPara$b_sigma)) 2 else priorPara$b_sigma

  Xbeta <- as.vector(X %*% beta)
  Zb <- as.vector(Z %*% bcoef)

  hazard_info <- cumulative_hazard_matrix(t = survObj$t, interval_bounds = priorPara$s,
                                          lambda = lambda)
  A <- hazard_info$A
  H0 <- hazard_info$H0
  h0 <- hazard_info$h0

  psd_b <- rep_len(proposal_sd_b, p_b)
  psd_beta <- rep_len(proposal_sd_beta, p_beta)
  psd_lambda <- rep_len(proposal_sd_lambda, J)
  psd_theta <- proposal_sd_theta

  save_iter <- seq.int(warmup + 1L, nIter, by = thin)
  nSave <- length(save_iter)
  save_map <- integer(nIter)
  save_map[save_iter] <- seq_len(nSave)

  beta.p <- matrix(NA_real_, nSave, p_beta)
  b.p <- matrix(NA_real_, nSave, p_b)
  lambda.p <- matrix(NA_real_, nSave, J)
  theta.p <- matrix(NA_real_, nSave, 1L)
  sigmaSq.p <- numeric(nSave);  sigmaStarSq.p <- numeric(nSave)
  etaSq.p <- numeric(nSave);    etaStarSq.p <- numeric(nSave)
  tauSq.p <- if (store_tau) matrix(NA_real_, nSave, p_b) else NULL
  tauStarSq.p <- if (store_tau) matrix(NA_real_, nSave, p_beta) else NULL
  loglik.p <- if (save_loglik == 1L) matrix(NA_real_, nSave, n) else NULL

  accept.b <- integer(p_b);  accept.beta <- integer(p_beta)
  accept.lambda <- integer(J);  accept.theta <- 0L

  pb <- NULL
  if (isTRUE(progress)) pb <- utils::txtProgressBar(min = 0, max = nIter, style = 3)

  for (M in seq_len(nIter)) {

    ## Periodic refresh of the incrementally updated quantities (Zb, Xbeta, H0, h0) from the current
    ## parameters, so that floating-point error cannot accumulate over long runs.
    if (refresh_every > 0L && M %% refresh_every == 0L) {
      Zb    <- as.vector(Z %*% bcoef)
      Xbeta <- as.vector(X %*% beta)
      hz    <- smcm_hazard_from_state(lambda, hazard_info)
      H0    <- hz$H0;  h0 <- hz$h0
    }

    ## ---- b ----------------------------------------------------------------##
    sd.b <- sqrt(sigmaSq * tauSq)
    if (length(unpen_b)) sd.b[unpen_b] <- intercept_sd
    H1   <- H0 * exp(Xbeta)
    logS <- smcfm_logS(H1, theta)
    ub <- UpdateRP.b.SMCFM.fast(Z, bcoef, Zb, logS, event, sd.b, psd_b)
    bcoef <- ub$b;  Zb <- ub$Zb

    ## ---- beta -------------------------------------------------------------##
    sd.beta <- sqrt(sigmaStarSq * tauStarSq)
    ube <- UpdateRP.beta.SMCFM.fast(X, beta, Xbeta, Zb, H0, theta, event, sd.beta, psd_beta)
    beta <- ube$beta;  Xbeta <- ube$Xbeta

    ## ---- lambda -----------------------------------------------------------##
    ul <- UpdateBH.lambda.SMCFM.fast(lambda, H0, h0, hazard_info, Xbeta, Zb, theta, event,
                                     priorPara$a, priorPara$b, psd_lambda)
    lambda <- ul$lambda;  H0 <- ul$H0;  h0 <- ul$h0

    ## ---- theta ------------------------------------------------------------##
    ut <- UpdateTheta.SMCFM.fast(theta, H0 * exp(Xbeta), Zb, event,
                                 priorPara$c, priorPara$d, psd_theta)
    theta <- ut$theta

    ## Acceptance rates are counted after warmup only
    if (M > warmup) {
      accept.b <- accept.b + ub$accept
      accept.beta <- accept.beta + ube$accept
      accept.lambda <- accept.lambda + ul$accept
      accept.theta <- accept.theta + ut$accept
    }

    ## Robbins-Monro adaptation during warmup only (frozen afterwards)
    if (isTRUE(adapt) && M <= warmup) {
      g <- 1 / (M + 10)^0.6
      psd_b      <- pmin(pmax(psd_b      * exp(g * (ub$accept  - target_accept)), 1e-3), 20)
      psd_beta   <- pmin(pmax(psd_beta   * exp(g * (ube$accept - target_accept)), 1e-3), 20)
      psd_lambda <- pmin(pmax(psd_lambda * exp(g * (ul$accept  - target_accept)), 1e-3), 20)
      psd_theta  <- min(max(psd_theta    * exp(g * (ut$accept  - target_accept)), 1e-3), 20)
    }

    ## ---- lasso hierarchy (shared with SMCM) -------------------------------##
    tauSq[pen_b] <- UpdateTau.fast(bcoef[pen_b], sigmaSq, etaSq)
    tauStarSq <- UpdateTauStar.fast(beta, sigmaStarSq, etaStarSq)
    sigmaSq <- UpdateSigma.fast(bcoef[pen_b], tauSq[pen_b], a_sigma, b_sigma)
    sigmaStarSq <- UpdateSigmaStar.fast(beta, tauStarSq, a_sigma, b_sigma)
    etaSq <- UpdateEta.fast(tauSq[pen_b], priorPara$r1, priorPara$delta1)
    etaStarSq <- UpdateEtaStar.fast(tauStarSq, priorPara$r2, priorPara$delta2)

    ## ---- store ------------------------------------------------------------==
    si <- save_map[M]
    if (si > 0L) {
      beta.p[si, ] <- beta;  b.p[si, ] <- bcoef
      lambda.p[si, ] <- lambda;  theta.p[si, 1L] <- theta
      sigmaSq.p[si] <- sigmaSq;  sigmaStarSq.p[si] <- sigmaStarSq
      etaSq.p[si] <- etaSq;  etaStarSq.p[si] <- etaStarSq
      if (store_tau) {
        tauSq.p[si, ] <- tauSq;  tauStarSq.p[si, ] <- tauStarSq
      }
      if (save_loglik == 1L) {
        loglik.p[si, ] <- as.numeric(logLik_SMCFM_fast(survObj, priorPara, bcoef, beta,
                                                       lambda, theta, hazard_info))
      }
    }

    if (!is.null(pb) && (M %% max(1L, floor(nIter / 10L)) == 0L || M == nIter)) {
      utils::setTxtProgressBar(pb, M)
    }
  }
  if (!is.null(pb)) close(pb)

  nPost <- nIter - warmup
  list(
    beta.p = beta.p, b.p = b.p, lambda.p = lambda.p, theta.p = theta.p,
    tauSq.p = tauSq.p, tauStarSq.p = tauStarSq.p,
    mcmcOutcome = list(
      initial = initial, priorPara = priorPara,
      sigmaSq.p = sigmaSq.p, sigmaStarSq.p = sigmaStarSq.p,
      etaSq.p = etaSq.p, etaStarSq.p = etaStarSq.p,
      accept.beta = accept.beta, accept.b = accept.b,
      accept.lambda = accept.lambda, accept.theta = accept.theta,
      acceptance.beta = accept.beta / nPost,            # post-warmup only
      acceptance.b = accept.b / nPost,
      acceptance.lambda = accept.lambda / nPost,
      acceptance.theta = accept.theta / nPost,
      proposal.sd = list(b = psd_b, beta = psd_beta, lambda = psd_lambda, theta = psd_theta),
      final = list(beta.ini = beta, b.ini = bcoef, lambda.ini = lambda, theta.ini = theta,
                   sigmaSq = sigmaSq, sigmaStarSq = sigmaStarSq,
                   etaSq = etaSq, etaStarSq = etaStarSq,
                   tauSq = tauSq, tauStarSq = tauStarSq,
                   Xbeta = Xbeta, Zb = Zb, H0 = H0, h0 = h0,
                   A = A, interval_id = hazard_info$interval_id)
    ),
    t = survObj$t,
    di = survObj$di,
    loglik.matrix = loglik.p,
    save_iter = save_iter,
    elapsed = NULL
  )
}
