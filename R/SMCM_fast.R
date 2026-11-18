# Optimized implementation for the Bayesian semiparametric
# mixture cure model (SMCM).
#
softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

# Piecewise-exponential baseline hazard using an interval-exposure matrix A.
#
# A[i, j] is the amount of time subject i spends in interval j.
# Hence H0 = A %*% lambda.
#
# Boundary convention matches the original cumulative_hazard():
# interval j is (interval_bounds[j], interval_bounds[j + 1]] for t > 0.
cumulative_hazard_matrix <- function(t, interval_bounds, lambda) {

  t <- as.numeric(t)
  interval_bounds <- as.numeric(interval_bounds)
  lambda <- as.numeric(lambda)

  n <- length(t)
  J <- length(lambda)

  if (length(interval_bounds) != J + 1L) {
    stop("length(interval_bounds) must equal length(lambda) + 1")
  }
  if (length(interval_bounds) < 2L || any(diff(interval_bounds) <= 0)) {
    stop("interval_bounds must be strictly increasing.")
  }
  if (any(!is.finite(t)) ||
      any(t < interval_bounds[1L]) ||
      any(t > tail(interval_bounds, 1L))) {
    stop("All observed times must lie inside the supplied interval boundaries.")
  }
  if (any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("All baseline hazard parameters lambda must be positive and finite.")
  }

  # Exposure matrix A: A[i, j] = time subject i spends in interval j
  A <- matrix(0, nrow = n, ncol = J)

  for (j in seq_len(J)) {
    A[, j] <- pmax(
      0,
      pmin(t, interval_bounds[j + 1L]) - interval_bounds[j]
    )
  }

  # Baseline cumulative hazard for every observation:
  # H0(t_i) = sum_j A[i,j] * lambda[j]
  H0 <- as.vector(A %*% lambda)

  # Find the interval containing each observed time.
  # This matches (s_j, s_{j+1}] used in the original code.
  interval_id <- findInterval(
    t,
    vec = interval_bounds,
    left.open = TRUE,
    rightmost.closed = TRUE
  )

  # Safety check
  interval_id <- pmin(pmax(interval_id, 1L), J)

  # Instantaneous piecewise-constant baseline hazard h0(t_i)
  h0 <- lambda[interval_id]

  list(
    h0 = h0,
    H0 = H0,
    A = A,
    interval_id = interval_id
  )
}

# Recompute H0 and h0 from a new lambda while reusing the fixed A matrix
# and interval memberships. This is useful outside the coordinate-wise
# lambda update (for example, validation or posterior likelihood evaluation).
smcm_hazard_from_state <- function(lambda, hazard_info) {
  H0 <- as.vector(hazard_info$A %*% lambda)
  h0 <- lambda[hazard_info$interval_id]
  list(H0 = H0, h0 = h0)
}

# b-dependent part of the observed-data log-likelihood.
smcm_loglik_b <- function(Zb, H1, delta) {
  event <- delta == 1
  cens <- !event
  out <- -sum(softplus(Zb))
  if (any(event)) out <- out + sum(Zb[event])
  if (any(cens)) out <- out + sum(softplus(Zb[cens] - H1[cens]))
  out
}

# beta-dependent part of the observed-data log-likelihood.
smcm_loglik_beta <- function(Xbeta, Zb, H0, delta) {
  H1 <- H0 * exp(Xbeta)
  event <- delta == 1
  cens <- !event
  out <- 0
  if (any(event)) out <- out + sum(Xbeta[event] - H1[event])
  if (any(cens)) out <- out + sum(softplus(Zb[cens] - H1[cens]))
  if (!is.finite(out)) -Inf else out
}

# lambda-dependent part of the observed-data log-likelihood.
smcm_loglik_lambda <- function(H0, h0, Xbeta, Zb, delta) {
  if (any(h0 <= 0) || any(!is.finite(h0))) return(-Inf)
  H1 <- H0 * exp(Xbeta)
  event <- delta == 1
  cens <- !event
  out <- 0
  if (any(event)) out <- out + sum(log(h0[event]) - H1[event])
  if (any(cens)) out <- out + sum(softplus(Zb[cens] - H1[cens]))
  if (!is.finite(out)) -Inf else out
}

# Full per-observation log-likelihood, useful for LOO/DIC and validation.
logLik_SMCM_fast <- function(survObj, priorPara, b.coef, beta.coef, lambda.coef,
                             hazard_info = NULL) {
  if (is.null(hazard_info)) {
    hazard_info <- cumulative_hazard_matrix(
      t = survObj$t,
      interval_bounds = priorPara$s,
      lambda = lambda.coef
    )
  }

  Zb<- as.vector(survObj$Z %*% b.coef)
  Xbeta <- as.vector(survObj$X %*% beta.coef)
  hz <- smcm_hazard_from_state(lambda.coef, hazard_info)
  H1 <- hz$H0 * exp(Xbeta)

  event <- survObj$di == 1
  cens <- !event
  ll <- -softplus(Zb)
  if (any(event)) {
    ll[event] <- ll[event] + Zb[event] + log(hz$h0[event]) +
      Xbeta[event] - H1[event]
  }
  if (any(cens)) {
    ll[cens] <- ll[cens] + softplus(Zb[cens] - H1[cens])
  }
  matrix(ll, nrow = 1L, dimnames = list(NULL, paste0("V", seq_along(ll))))
}

UpdateRP.b.fast <- function(Z, b, Zb, H1, delta, sd.b,
                            proposal.sd = rep(1, length(b)) ){
  p <- length(b)
  if (length(proposal.sd) == 1L) proposal.sd <- rep(proposal.sd, p)
  accept <- integer(p)
  current.loglik <- smcm_loglik_b(Zb, H1, delta)

  for (k in seq_len(p)) {
    b.prop <- rnorm(1L, mean = b[k], sd = proposal.sd[k])
    diff <- b.prop - b[k]
    Zb.prop <- Zb + Z[, k] * diff
    prop.loglik <- smcm_loglik_b(Zb.prop, H1, delta)

    logR <- prop.loglik - current.loglik +
      dnorm(b.prop, 0, sd.b[k], log = TRUE) -
      dnorm(b[k], 0, sd.b[k], log = TRUE)

    if (!is.na(logR) && log(runif(1L)) < logR) {   
      b[k] <- b.prop
      Zb <- Zb.prop
      current.loglik <- prop.loglik
      accept[k] <- 1L
    }
  }

  list(b = b, Zb = Zb, loglik = current.loglik, accept = accept)
}

UpdateRP.beta.fast <- function(X, beta, Xbeta, Zb, H0, delta, sd.beta,
                               proposal.sd = rep(1, length(beta))) {
  p <- length(beta)
  if (length(proposal.sd) == 1L) proposal.sd <- rep(proposal.sd, p)
  accept <- integer(p)
  current.loglik <- smcm_loglik_beta(Xbeta, Zb, H0, delta)
  H1 <- H0 * exp(Xbeta)

  for (k in seq_len(p)) {
    beta.prop <- rnorm(1L, mean = beta[k], sd = proposal.sd[k])
    diff <- beta.prop - beta[k]
    Xbeta.prop <- Xbeta + X[, k] * diff
    H1.prop <- H0 * exp(Xbeta.prop)

    event <- delta == 1
    cens <- !event
    prop.loglik <- 0
    if (any(event)) prop.loglik <- prop.loglik + sum(Xbeta.prop[event] - H1.prop[event])
    if (any(cens)) prop.loglik <- prop.loglik + sum(softplus(Zb[cens] - H1.prop[cens]))
    if (!is.finite(prop.loglik)) prop.loglik <- -Inf

    logR <- prop.loglik - current.loglik +
      dnorm(beta.prop, 0, sd.beta[k], log = TRUE) -
      dnorm(beta[k], 0, sd.beta[k], log = TRUE)

    if (!is.na(logR) && log(runif(1L)) < logR) {
      beta[k] <- beta.prop
      Xbeta <- Xbeta.prop
      H1 <- H1.prop
      current.loglik <- prop.loglik
      accept[k] <- 1L
    }
  }

  list(beta = beta, Xbeta = Xbeta, H1 = H1,
       loglik = current.loglik, accept = accept)
}
#'
#'
UpdateBH.lambda.fast <- function(lambda, H0, h0, hazard_info, Xbeta, Zb, delta,
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
  current.loglik <- smcm_loglik_lambda(H0, h0, Xbeta, Zb, delta)
  
  
  for (k in seq_len(J)) {
    
    lambda.current <- lambda[k]
    ## ----------------------------------------------------------
    ## 1. Log-scale Gaussian random-walk proposal
    ##
    ## log(lambda_k*) =
    ##   log(lambda_k) + epsilon
    ##
    ## epsilon ~ N(0, proposal.sd[k]^2)
    ## ----------------------------------------------------------
    log.lambda.current <- log(lambda.current)
    log.lambda.prop <- log.lambda.current + rnorm(1L, mean = 0, sd = proposal.sd[k])
    
    ## avoid overflow before exponentiation
    if (!is.finite(log.lambda.prop)) {
      next
    }
    lambda.prop <- exp(log.lambda.prop)
    if (!is.finite(lambda.prop) || lambda.prop <= 0) {
      next
    }
    
    ## ----------------------------------------------------------
    ## 2. Incremental update of cumulative baseline hazard
    ##
    ## H0* = H0 + A[,k] (lambda_k* - lambda_k)
    ## ----------------------------------------------------------
    
    diff <- lambda.prop - lambda.current
    H0.prop <- H0 + hazard_info$A[, k] * diff
    if (any(!is.finite(H0.prop)) || any(H0.prop < 0)) {
      next
    }
    ## ----------------------------------------------------------
    ## 3. Instantaneous baseline hazard
    ##
    ## Only observations belonging to interval k change.
    ## ----------------------------------------------------------
    h0.prop <- h0
    idx <- hazard_info$interval_id == k
    if (any(idx)) { h0.prop[idx] <- lambda.prop  }
    
    ## ----------------------------------------------------------
    ## 4. Proposed lambda-dependent log-likelihood
    ## ----------------------------------------------------------
    prop.loglik <- smcm_loglik_lambda(H0.prop, h0.prop, Xbeta, Zb, delta)
    if (!is.finite(prop.loglik)) {
      next
    }
    
    ## ----------------------------------------------------------
    ## 5. Gamma prior for lambda_k
    ## lambda_k ~ Gamma(a, b)
    ## b is the RATE parameter.
    ## ----------------------------------------------------------
    logprior.current <- dgamma(lambda.current, shape = a, rate = b, log = TRUE )
    logprior.prop <- dgamma( lambda.prop, shape = a, rate = b, log = TRUE)
    
    ## ----------------------------------------------------------
    ## 6. Proposal-density correction
    ##
    ## Proposal:
    ##
    ## log(lambda*) | log(lambda)
    ##     ~ N(log(lambda), proposal.sd^2)
    ##
    ## On the lambda scale this is log-normal.
    ##
    ## q(lambda | lambda*) / q(lambda* | lambda)
    ##     = lambda* / lambda
    ##
    ## Therefore:
    ##
    ## log q(lambda | lambda*) -
    ## log q(lambda* | lambda)
    ##
    ## = log(lambda*) - log(lambda)
    ## ----------------------------------------------------------
    
    
    ## ----------------------------------------------------------
    ## 7. Metropolis-Hastings log-ratio
    ## ----------------------------------------------------------
    
    logR <- prop.loglik - current.loglik + logprior.prop - logprior.current + 
      log.lambda.prop - log.lambda.current
    
    ## ----------------------------------------------------------
    ## 8. Accept/reject
    ## ----------------------------------------------------------
    if (!is.na(logR) && log(runif(1L)) < logR) {
      lambda[k] <- lambda.prop
      H0 <- H0.prop
      h0 <- h0.prop
      current.loglik <- prop.loglik  ## Carry current likelihood forward.
      accept[k] <- 1L
    }
    
  }

  list(lambda = lambda, H0 = H0, h0 = h0, loglik = current.loglik, accept = accept)
  
}

# Bayesian-Lasso hierarchy. A small floor avoids overflow in 1 / coef^2 when
# a coefficient is numerically zero; it does not alter ordinary nonzero values.
UpdateTau.fast <- function(b, sigmaSq, etaSq, eps = 1e-12) {
  absb <- pmax(abs(b), eps)
  nu <- sqrt(etaSq * sigmaSq) / absb
  gam <- vapply(nu, function(mu) {
    g <- SuppDists::rinvGauss(1L, nu = mu, lambda = etaSq)
    while (!is.finite(g) || g <= 0) {
      g <- SuppDists::rinvGauss(1L, nu = mu, lambda = etaSq)
    }
    g
  }, numeric(1L))
  1 / gam
}

UpdateTauStar.fast <- function(beta, sigmaStarSq, etaStarSq, eps = 1e-12) {
  absb <- pmax(abs(beta), eps)
  nu <- sqrt(etaStarSq * sigmaStarSq) / absb
  gam <- vapply(nu, function(mu) {
    g <- SuppDists::rinvGauss(1L, nu = mu, lambda = etaStarSq)
    while (!is.finite(g) || g <= 0) {
      g <- SuppDists::rinvGauss(1L, nu = mu, lambda = etaStarSq)
    }
    g
  }, numeric(1L))
  1 / gam
}

UpdateSigma.fast <- function(b, tauSq, a_sigma = 3, b_sigma = 2) {
  shape <- a_sigma + length(b) / 2
  rate <- b_sigma + sum((b^2) / (2 * tauSq))
  1 / rgamma(1L, shape = shape, rate = rate)      # inverse-gamma(shape, scale = rate)
}

UpdateSigmaStar.fast <- function(beta, tauStarSq,  a_sigma = 3, b_sigma = 2) {
  shape <- a_sigma + length(beta) / 2
  rate <- b_sigma + sum((beta^2) / (2 * tauStarSq))
  1 / rgamma(1L, shape = shape, rate = rate)      # inverse-gamma(shape, scale = rate)
}

UpdateEta.fast <- function(tauSq, r1, delta1) {
  rgamma(1L, shape = r1 + length(tauSq), rate = delta1 + sum(tauSq) / 2)
}

UpdateEtaStar.fast <- function(tauStarSq, r2, delta2) {
  rgamma(1L, shape = r2 + length(tauStarSq), rate = delta2 + sum(tauStarSq) / 2)
}

#' Optimized single-chain R reference implementation.
#' @title Function to MCMC algorithm for Semi-parametric Mixture Cure Model (SMCM)
#'
#' @description
#' This function obtains posterior samples for all model parameters using an MCMC algorithm 
#' for a single chain when applying the \code{MCMC.SMCM.fast} framework in 
#' Kızılaslan and Vitelli (2025) to the data.
#'
#' @name MCMC.SMCM.fast
#' @export
MCMC.SMCM.fast <- function(survObj, priorPara, initial,
                           nIter = 6000L, warmup = 1000L, thin = 1L,
                           save_loglik = 0L, store_tau = TRUE,
                           proposal_sd_b = 1,
                           proposal_sd_beta = 1,
                           proposal_sd_lambda = 0.5,
                           adapt = TRUE, target_accept = 0.35,
                           progress = interactive()) {
  stopifnot(nIter >= 1L, warmup >= 0L, warmup < nIter, thin >= 1L)

  X <- as.matrix(survObj$X)
  Z <- as.matrix(survObj$Z)
  delta <- as.integer(survObj$di)
  n <- nrow(X)
  p_beta <- ncol(X)
  p_b <- ncol(Z)
  J <- priorPara$J

  if (length(priorPara$s) - 1L != J) {
    stop("priorPara$J does not match length(priorPara$s)-1.")
  }

  beta <- as.numeric(initial$beta.ini)
  bcoef <- as.numeric(initial$b.ini)
  lambda <- as.numeric(initial$lambda.ini)
  sigmaSq <- initial$sigmaSq
  sigmaStarSq <- initial$sigmaStarSq
  etaSq <- initial$etaSq
  etaStarSq <- initial$etaStarSq
  tauSq <- as.numeric(initial$tauSq)

  ## Intercept handling: if Z[, 1] is an intercept column and penalize_intercept
  ## is FALSE, it gets a vague N(0, intercept_sd^2) prior and is excluded from the
  ## lasso hierarchy (tau, sigma, eta updates use the slopes only).
  has_int      <- p_b > 0L && all(Z[, 1L] == 1)
  unpen_b      <- if (!isTRUE(priorPara$penalize_intercept) && has_int) 1L else integer(0)
  pen_b        <- setdiff(seq_len(p_b), unpen_b)
  intercept_sd <- if (is.null(priorPara$intercept_sd)) 10 else priorPara$intercept_sd
  tauSq[unpen_b] <- 1

  ## Inverse-gamma hyperparameters of sigma^2 (previously hard-coded 3 and 2)
  a_sigma <- if (is.null(priorPara$a_sigma)) 3 else priorPara$a_sigma
  b_sigma <- if (is.null(priorPara$b_sigma)) 2 else priorPara$b_sigma
  tauStarSq <- as.numeric(initial$tauStarSq)

  Xbeta <- as.vector(X %*% beta)
  Zb <- as.vector(Z %*% bcoef)

  # Build the interval-exposure matrix A once.  A and interval_id remain
  # fixed throughout the MCMC because the observed times and cutpoints are fixed.
  hazard_info <- cumulative_hazard_matrix(
    t = survObj$t,
    interval_bounds = priorPara$s,
    lambda = lambda
  )
  A <- hazard_info$A
  H0 <- hazard_info$H0
  h0 <- hazard_info$h0

  psd_b <- rep_len(proposal_sd_b, p_b)
  psd_beta <- rep_len(proposal_sd_beta, p_beta)
  psd_lambda <- rep_len(proposal_sd_lambda, J)

  save_iter <- seq.int(warmup + 1L, nIter, by = thin)
  nSave <- length(save_iter)
  save_map <- integer(nIter)
  save_map[save_iter] <- seq_len(nSave)

  beta.p <- matrix(NA_real_, nSave, p_beta)
  b.p <- matrix(NA_real_, nSave, p_b)
  lambda.p <- matrix(NA_real_, nSave, J)
  sigmaSq.p <- numeric(nSave)
  sigmaStarSq.p <- numeric(nSave)
  etaSq.p <- numeric(nSave)
  etaStarSq.p <- numeric(nSave)
  tauSq.p <- if (store_tau) matrix(NA_real_, nSave, p_b) else NULL
  tauStarSq.p <- if (store_tau) matrix(NA_real_, nSave, p_beta) else NULL
  loglik.p <- if (save_loglik == 1L) matrix(NA_real_, nSave, n) else NULL

  accept.b <- integer(p_b)
  accept.beta <- integer(p_beta)
  accept.lambda <- integer(J)

  pb <- NULL
  if (isTRUE(progress)) pb <- utils::txtProgressBar(min = 0, max = nIter, style = 3)

  for (M in seq_len(nIter)) {
    sd.b <- sqrt(sigmaSq * tauSq)
    if (length(unpen_b)) sd.b[unpen_b] <- intercept_sd      # vague prior on the intercept
    H1 <- H0 * exp(Xbeta)
    ub <- UpdateRP.b.fast(Z, bcoef, Zb, H1, delta, sd.b, psd_b)
    bcoef <- ub$b
    Zb <- ub$Zb

    sd.beta <- sqrt(sigmaStarSq * tauStarSq)
    ube <- UpdateRP.beta.fast(X, beta, Xbeta, Zb, H0, delta,
                              sd.beta, psd_beta)
    beta <- ube$beta
    Xbeta <- ube$Xbeta

    ul <- UpdateBH.lambda.fast(lambda, H0, h0, hazard_info, Xbeta, Zb, delta,
                               priorPara$a, priorPara$b, psd_lambda)
    lambda <- ul$lambda
    H0 <- ul$H0
    h0 <- ul$h0

    ## Acceptance rates are counted after warmup only
    if (M > warmup) {
      accept.b      <- accept.b + ub$accept
      accept.beta   <- accept.beta + ube$accept
      accept.lambda <- accept.lambda + ul$accept
    }

    ## Robbins-Monro adaptation of the coordinate-wise proposal SDs during
    ## warmup only (frozen afterwards, so the saved chain targets the posterior).
    if (isTRUE(adapt) && M <= warmup) {
      g <- 1 / (M + 10)^0.6
      psd_b      <- pmin(pmax(psd_b      * exp(g * (ub$accept  - target_accept)), 1e-3), 20)
      psd_beta   <- pmin(pmax(psd_beta   * exp(g * (ube$accept - target_accept)), 1e-3), 20)
      psd_lambda <- pmin(pmax(psd_lambda * exp(g * (ul$accept  - target_accept)), 1e-3), 20)
    }

    tauSq[pen_b] <- UpdateTau.fast(bcoef[pen_b], sigmaSq, etaSq)
    tauStarSq <- UpdateTauStar.fast(beta, sigmaStarSq, etaStarSq)
    sigmaSq <- UpdateSigma.fast(bcoef[pen_b], tauSq[pen_b], a_sigma, b_sigma)
    sigmaStarSq <- UpdateSigmaStar.fast(beta, tauStarSq, a_sigma, b_sigma)
    etaSq <- UpdateEta.fast(tauSq[pen_b], priorPara$r1, priorPara$delta1)
    etaStarSq <- UpdateEtaStar.fast(tauStarSq, priorPara$r2, priorPara$delta2)

    si <- save_map[M]
    if (si > 0L) {
      beta.p[si, ] <- beta
      b.p[si, ] <- bcoef
      lambda.p[si, ] <- lambda
      sigmaSq.p[si] <- sigmaSq
      sigmaStarSq.p[si] <- sigmaStarSq
      etaSq.p[si] <- etaSq
      etaStarSq.p[si] <- etaStarSq
      if (store_tau) {
        tauSq.p[si, ] <- tauSq
        tauStarSq.p[si, ] <- tauStarSq
      }
      if (save_loglik == 1L) {
        loglik.p[si, ] <- as.numeric(logLik_SMCM_fast(
          survObj, priorPara, bcoef, beta, lambda, hazard_info
        ))
      }
    }

    if (!is.null(pb) && (M %% max(1L, floor(nIter / 10L)) == 0L || M == nIter)) {
      utils::setTxtProgressBar(pb, M)
    }
  }
  if (!is.null(pb)) close(pb)

  list(
    beta.p = beta.p,
    b.p = b.p,
    lambda.p = lambda.p,
    tauSq.p = tauSq.p,
    tauStarSq.p = tauStarSq.p,
    mcmcOutcome = list(
      initial = initial,
      priorPara = priorPara,
      sigmaSq.p = sigmaSq.p,
      sigmaStarSq.p = sigmaStarSq.p,
      etaSq.p = etaSq.p,
      etaStarSq.p = etaStarSq.p,
      accept.beta = accept.beta,
      accept.b = accept.b,
      accept.lambda = accept.lambda,
      acceptance.beta = accept.beta / (nIter - warmup),     # post-warmup only
      acceptance.b = accept.b / (nIter - warmup),
      acceptance.lambda = accept.lambda / (nIter - warmup),
      proposal.sd = list(b = psd_b, beta = psd_beta, lambda = psd_lambda),
      final = list(beta.ini = beta, b.ini = bcoef, lambda.ini = lambda,
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
#'