#' @title  Fit Bayesian Hierarchical Semi-parametric Mixture Cure Frailty Model (HSMCFM) with RStan
#'
#' @description
#' Fits the Bayesian HSMCFM model with \pkg{RStan}, using hierarchical (Bayesian-lasso type) priors
#' for \eqn{\bm b} and \eqn{\bm \beta} as in Kızılaslan and Vitelli (2025). The uncured subjects have a
#' gamma frailty (mean 1, variance 1/theta).
#'
#' @name fit.HSMCFM.RStan
#' @import rstan
#'
#' @param data list with \code{observed_time}, \code{delta}, \code{X} (latency part, n x p1, coefficients beta,
#' NO intercept) and \code{Z} (incidence/cure part, n x p2, coefficients b, FIRST column = intercept).
#' @param hyperpar list with ALL of \code{r1, delta1, r2, delta2, a, b, c, d, a_sigma, b_sigma, a_sigmaStar, b_sigmaStar}.
#' Gamma priors are (shape, rate): eta^2 (r1, delta1), eta*^2 (r2, delta2), lambda (a, b), theta (c, d).
#' Inverse-gamma priors are (shape, scale): sigma^2 (a_sigma, b_sigma), sigma*^2 (a_sigmaStar, b_sigmaStar).
#' @param nchains,nIter,warmup,thin MCMC settings (nIter includes warmup).
#' @param standardize if TRUE, X and Z are centred/scaled (see \code{scale.info}, \code{scale.apply}); the model is
#' fitted on the scaled data and the draws in \code{b.chains}, \code{beta.chains}, \code{lambda.chains} are
#' transformed back to the ORIGINAL covariate scale (theta is not affected by the scaling).
#' @param scale.binary if TRUE, binary dummy columns are scaled too (the intercept never is).
#' @param probs quantiles of the event times used as inner cut points (0 = single interval).
#' @param save_loglik 0/1, save the observation-level log-likelihood (needed for DIC/LPML/LOO).
#' @param seed random seed.
#' @param adapt_delta,max_treedepth Stan control. If diagnostics fail the model is refitted once
#' with \code{adapt_delta = max(adapt_delta, 0.95)}.
#' @param rhat_threshold Rhat threshold used in \code{check.stan.rhat}.
#' @param stan.file optional path to the Stan file (default: HSMCFM.stan of the package).
#' @param refresh Stan progress interval.
#' @param penalize_intercept FALSE (default, as in the Rcpp samplers): Z[, 1] is an intercept with a N(0, intercept_sd^2)
#' prior, outside the lasso hierarchy. TRUE: the intercept is penalised like the other coefficients.
#' @param intercept_sd prior sd of the unpenalised intercept.
#' @param noncentered TRUE/FALSE: non-centred parameterisation of b and beta in Stan (better geometry when p is large
#' or the data are weakly informative). NULL (default) = automatic: TRUE if p1 + p2 >= 20 (heuristic, check with a pilot run).
#'
#' @return list with \code{fit} (stanfit, fitted/scaled scale), \code{fit.results} (posterior means, ORIGINAL scale),
#' \code{b.chains}, \code{beta.chains}, \code{lambda.chains} (ORIGINAL scale), \code{theta.chains},
#' \code{scaled.chains}, \code{other.chains}, \code{Rhat}, \code{diagnostics}, \code{loo.results},
#' \code{loglik.matrix.chains} (list of draws x n matrices, one per chain), \code{data} (as FITTED),
#' \code{data.orig}, \code{scale.info}, \code{priors}, \code{mcmc.info}.
#'
#' @examples
#' n <- 300
#' b.true <- c( 0.4, 0.5, 0.1)
#' beta.true <- c(1, 0.2 )
#' baseline.hazard.rates <- 1
#' intervals <- c(0, 16)
#' dat1 <- simSMCM(n, b.true, beta.true, baseline.hazard.rates, intervals, seed = 2026,
#'   cens.start = 30/365, cens.end = 32,  prob.cov.X = 0.5, prob.cov.Z = 0.5, same.cov = TRUE )
#' priorPar = list( r1 = 1, delta1 = 1, r2 = 1, delta2 = 1, a = 1, b = 1,
#'                  c = 1, d = 1,
#'                  a_sigma = 3, b_sigma = 2, a_sigmaStar = 3, b_sigmaStar = 2 )
#' out <- fit.HSMCFM.RStan(data = dat1, hyperpar = priorPar, nchains = 2, nIter = 2500,
#'   warmup = 500, thin = 1, standardize = FALSE, probs = 0, save_loglik = 1, seed = 2026)
#'
#' @references
#' Kızılaslan, F., Vitelli, V. (2025). Bayesian Semiparametric Mixture Cure (Frailty) Model, arXiv \url{https://arxiv.org/abs/2512.08173}.
#'
#' @export
fit.HSMCFM.RStan <- function(data, hyperpar = list(), nchains = 2, nIter = 6000, warmup = 1000, thin = 1,
                             standardize = TRUE, scale.binary = TRUE, probs = 0, save_loglik = 1, seed = NULL,
                             adapt_delta = 0.80, max_treedepth = 10, rhat_threshold = 1.1,
                             stan.file = NULL, refresh = NULL, noncentered = NULL,
                             penalize_intercept = FALSE, intercept_sd = 10) {

  ## ---------------------------------------------------------------- checks
  required <- c("r1", "delta1", "r2", "delta2", "a", "b", "c", "d",
                "a_sigma", "b_sigma", "a_sigmaStar", "b_sigmaStar")
  if (is.null(hyperpar$a_sigma))     hyperpar$a_sigma     <- 3                     # defaults as in the Rcpp samplers
  if (is.null(hyperpar$b_sigma))     hyperpar$b_sigma     <- 2
  if (is.null(hyperpar$a_sigmaStar)) hyperpar$a_sigmaStar <- hyperpar$a_sigma    # same prior for sigma* by default
  if (is.null(hyperpar$b_sigmaStar)) hyperpar$b_sigmaStar <- hyperpar$b_sigma
  miss <- setdiff(required, names(hyperpar))
  if (length(miss) > 0) stop("Missing hyperparameters: ", paste(miss, collapse = ", "))
  if (any(!is.finite(unlist(hyperpar[required]))) || any(unlist(hyperpar[required]) <= 0))
    stop("All hyperparameters must be positive and finite.")
  if (warmup >= nIter) stop("warmup must be smaller than nIter.")
  save_loglik <- as.integer(save_loglik)

  t     <- as.numeric(data$observed_time)
  delta <- as.numeric(data$delta)
  n     <- length(t)
  if (anyNA(t) || anyNA(delta) || length(delta) != n) stop("observed_time/delta: NA or length mismatch.")
  if (any(t < 0)) stop("observed_time must be non-negative.")
  if (!all(delta %in% c(0, 1))) stop("delta must be 0/1.")

  X0 <- as.matrix(data$X)
  Z0 <- as.matrix(data$Z)
  if (nrow(X0) != n || nrow(Z0) != n) stop("X and Z must have n rows.")
  if (anyNA(X0) || anyNA(Z0)) stop("X and Z contain NA.")
  p1 <- ncol(X0); p2 <- ncol(Z0)
  if (is.null(noncentered)) noncentered <- (p1 + p2 >= 20)
  noncentered <- as.integer(isTRUE(as.logical(noncentered)))
  if ((isTRUE(standardize) || !isTRUE(penalize_intercept)) && !all(Z0[, 1] == 1))
    stop("The first column of Z must be the intercept (a column of ones), or use penalize_intercept = TRUE (and standardize = FALSE).")

  if (!is.null(seed)) set.seed(seed)
  seed.sampling <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else as.integer(seed)

  ## ------------------------------------------------------------ covariates
  info.X <- info.Z <- NULL
  if (isTRUE(standardize)) {
    info.X <- scale.info(X0, scale.binary);  X <- scale.apply(X0, info.X)
    info.Z <- scale.info(Z0, scale.binary);  Z <- scale.apply(Z0, info.Z)
  } else {
    X <- X0
    Z <- Z0
  }
  survObj <- list(t = t, di = delta, X = X, Z = Z)          # data AS FITTED (scaled if standardize = TRUE)

  ## ------------------------------------------------------------ cut points
  ## Intervals (s_j, s_{j+1}]; last cut point = max(t) (no exposure beyond max(t), posterior unchanged).
  event.time.points <- sort(t[delta == 1])
  inner <- if (is.null(probs) || (length(probs) == 1 && probs == 0)) {
    numeric(0)
  } else {
    if (any(probs <= 0 | probs >= 1)) stop("probs must lie in (0, 1) (or be 0 for J = 1).")
    as.numeric(quantile(event.time.points, probs = sort(probs), names = FALSE))
  }
  s  <- c(0, inner, max(t))
  if (any(diff(s) <= 0))
    stop("Cut points are not strictly increasing (ties in event times?). Change 'probs'.")
  JJ <- length(s) - 1

  ## ------------------------------------------------------------- priors
  priorPara <- list(s = s, J = JJ)
  for (nm in required) priorPara[[nm]] <- hyperpar[[nm]]

  data.stan.model <- list(
    n = n, p1 = p1, p2 = p2, JJ = JJ, s = s, t = t, delta = delta, X = X, Z = Z,
    lambda_a0    = priorPara$a,           lambda_b0    = priorPara$b,
    theta_a0     = priorPara$c,           theta_b0     = priorPara$d,
    sigma_a0     = priorPara$a_sigma,     sigma_b0     = priorPara$b_sigma,
    sigmaStar_a0 = priorPara$a_sigmaStar, sigmaStar_b0 = priorPara$b_sigmaStar,
    r1 = priorPara$r1, delta1 = priorPara$delta1,
    r2 = priorPara$r2, delta2 = priorPara$delta2,
    save_loglik = save_loglik, noncentered = noncentered,
    penalize_intercept = as.integer(isTRUE(penalize_intercept)), intercept_sd = intercept_sd
  )

  ## ------------------------------------------------------------ sampling
  HSMCFM.model <- .get.hsmcfm.model(stan.file)              # compiled once per session
  if (is.null(refresh)) refresh <- max(1L, floor(nIter / 2))
  n.cores <- max(1L, min(nchains, parallel::detectCores(), na.rm = TRUE))

  sampler.args <- list(object = HSMCFM.model, data = data.stan.model,
                       iter = nIter, warmup = warmup, thin = thin, chains = nchains,
                       cores = n.cores, refresh = refresh, seed = seed.sampling)
  if (save_loglik != 1) {
    sampler.args$pars    <- "log_lik"
    sampler.args$include <- FALSE
  }
  run.sampler <- function(ad)
    do.call(rstan::sampling,
            c(sampler.args, list(control = list(adapt_delta = ad, stepsize = 0.1,
                                                max_treedepth = max_treedepth))))

  pars.main <- c("b", "beta", "lambda", "theta", "sigmaSq", "etaSq", "sigmaStarSq", "etaStarSq")
  pars.diag <- c(pars.main, "lp__")
  get.diag <- function(fit) {
    sm <- rstan::summary(fit, pars = pars.diag)$summary
    list(summary  = sm,
         rhat.ok  = check.stan.rhat(fit, rhat_threshold = rhat_threshold, pars = pars.diag),
         max.rhat = max(sm[, "Rhat"], na.rm = TRUE),
         min.neff = min(sm[, "n_eff"], na.rm = TRUE),
         n.div    = rstan::get_num_divergent(fit))
  }
  is.ok <- function(d) d$rhat.ok && d$n.div == 0

  fit     <- run.sampler(adapt_delta)
  diag    <- get.diag(fit)
  elapsed <- sum(rstan::get_elapsed_time(fit))
  refit   <- FALSE
  if (!is.ok(diag)) {
    fit     <- run.sampler(max(adapt_delta, 0.95))          # HIGHER adapt_delta on the retry
    diag    <- get.diag(fit)
    elapsed <- elapsed + sum(rstan::get_elapsed_time(fit))
    refit   <- TRUE
    if (!is.ok(diag))
      warning(sprintf("Diagnostics still not clean after refit: max Rhat = %.3f, divergences = %d.",
                      diag$max.rhat, diag$n.div))
  }

  ## ------------------------------------------------- posterior (fitted scale)
  post      <- rstan::extract(fit, pars = pars.main)
  b.sc      <- as.matrix(post$b)
  beta.sc   <- as.matrix(post$beta)
  lambda.sc <- as.matrix(post$lambda)
  theta.mcmc.chains <- as.numeric(post$theta)

  ## ---------------------------------- DIC / LPML / LOO: ALWAYS on the fitted scale
  loo.results <- NULL
  loglik.matrix.chains <- NULL
  if (save_loglik == 1) {
    ll.arr <- rstan::extract(fit, pars = "log_lik", permuted = FALSE)    # iterations x chains x n
    loglik.matrix.chains <- lapply(seq_len(dim(ll.arr)[2]),
                                   function(ch) matrix(ll.arr[, ch, ], nrow = dim(ll.arr)[1]))
    loo.results <- mcmc.dic.lpml(survObj, priorPara, loglik.matrix.chains,
                                 b.hat.mcmc      = colMeans(b.sc),
                                 beta.hat.mcmc   = colMeans(beta.sc),
                                 lambda.hat.mcmc = colMeans(lambda.sc),
                                 stan.model = TRUE, frailty = TRUE,
                                 theta.hat.mcmc  = mean(theta.mcmc.chains))
  }

  ## ------------------------------------ back to the ORIGINAL covariate scale
  ## backtransform.b/beta/lambda are the helpers from the HSMCM code (theta needs no transformation:
  ## the frailty is multiplicative on the hazard and independent of the covariate scaling).
  if (isTRUE(standardize)) {
    b.orig      <- backtransform.b(b.sc, info.Z)
    beta.orig   <- backtransform.beta(beta.sc, info.X)
    lambda.orig <- backtransform.lambda(lambda.sc, beta.orig, info.X)
  } else {
    b.orig <- b.sc; beta.orig <- beta.sc; lambda.orig <- lambda.sc
  }

  pz.hat <- mean(plogis(Z %*% t(b.sc)))      # invariant to the scaling of Z

  fit.results <- c(b.hat       = colMeans(b.orig),
                   beta.hat    = colMeans(beta.orig),
                   lambda.hat  = colMeans(lambda.orig),
                   theta.hat   = mean(theta.mcmc.chains),
                   sigmaSq.hat = mean(post$sigmaSq),                     # fitted scale
                   sigmaStarSq.hat = mean(post$sigmaStarSq),             # fitted scale
                   pz.hat      = pz.hat,
                   time        = elapsed / 60)

  ret <- list(fit = fit,
              fit.results = fit.results,
              b.chains = b.orig,
              beta.chains = beta.orig,
              lambda.chains = lambda.orig,
              theta.chains = theta.mcmc.chains,
              scaled.chains = list(b = b.sc, beta = beta.sc, lambda = lambda.sc),
              other.chains = post[c("sigmaSq", "etaSq", "sigmaStarSq", "etaStarSq")],
              Rhat = diag$summary[, "Rhat"],
              diagnostics = list(max.rhat = diag$max.rhat, min.neff = diag$min.neff,
                                 n.divergent = diag$n.div, refitted = refit),
              loo.results = loo.results,
              loglik.matrix.chains = loglik.matrix.chains,
              data = survObj,
              data.orig = list(t = t, di = delta, X = X0, Z = Z0),
              scale.info = list(X = info.X, Z = info.Z),
              priors = priorPara,
              mcmc.info = list(nchains = nchains, nIter = nIter, warmup = warmup, thin = thin,
                                noncentered = as.logical(noncentered),
                                penalize_intercept = isTRUE(penalize_intercept), intercept_sd = intercept_sd))
  class(ret) <- "HSMCFM.RStan"
  ret
}

## ------------------------------------------------------------------------------------------
## Compile the Stan model once per session (and again only if the .stan file changes).
.hsmcfm.cache <- new.env(parent = emptyenv())

#' @noRd
.get.hsmcfm.model <- function(stan.file = NULL) {
  if (is.null(stan.file)) {
    stan.file <- system.file("stan", "HSMCFM.stan", package = "BayesSMCM")
    if (!nzchar(stan.file)) stan.file <- "~/BayesSMCM/inst/stan/HSMCFM.stan"
  }
  stan.file <- normalizePath(stan.file, mustWork = TRUE)
  key <- paste(stan.file, as.character(file.mtime(stan.file)))
  if (is.null(.hsmcfm.cache[[key]]))
    .hsmcfm.cache[[key]] <- rstan::stan_model(file = stan.file)
  .hsmcfm.cache[[key]]
}
