#' fit.MCMC.SMCM
#' @title Fit Bayesian Semi-parametric Mixture Cure Model (SMCM) via MCMC algorithm
#'
#' @description
#' This is the main function to fit the Bayesian SMCM model using MCMC algorithm
#' in Kızılaslan and Vitelli (2025).
#' Fits the Bayesian semiparametric mixture cure model designed to work with
#' MCMC.SMCM.fast.
#'
#' @name fit.MCMC.SMCFM
#'
#' @param data A list containing observed_time, delta, X and Z.
#' @param hyperpar Prior parameters r1, delta1, r2, delta2, a and b.
#' @param nchains Number of MCMC chains.
#' @param nIter Total number of MCMC iterations per chain.
#' @param warmup Number of warmup iterations. Warmup is handled inside
#'   MCMC.SMCM.fast.
#' @param thin Storage/thinning interval after warmup. Thinning is handled
#'   inside MCMC.SMCM.fast.
#' @param mcmc.parallel One of "lapply", "parLapply" or "foreach".
#' @param standardize Logical. If \code{TRUE} (default), the columns of X and Z
#'   are centred and scaled to unit standard deviation before the MCMC, as in
#'   glmnet; the intercept column of Z (any constant column) is never scaled.
#'   If \code{FALSE}, X and Z are used unchanged. See Note.
#' @param scale.binary Logical. Only used when \code{standardize = TRUE}. If
#'   \code{TRUE} (default), binary (integer-valued) columns are centred and
#'   scaled like the continuous ones; if \code{FALSE}, they are left unchanged.
#' @param probs Quantile probabilities used for internal piecewise-exponential
#'   cut points; use 0 for a single interval.
#' @param save_loglik Integer/logical indicating whether observation-level
#'   log-likelihood values are saved.
#' @param seed Optional random seed.
#' @param store_tau Logical; whether local shrinkage scales are stored.
#' @param proposal_sd_b Gaussian random-walk proposal SD for b.
#' @param proposal_sd_beta Gaussian random-walk proposal SD for beta.
#' @param proposal_sd_lambda Log-Gaussian random-walk proposal SD for lambda.
#' @param init_method Initial coefficient strategy: "legacy" reproduces the
#'   original chain-specific 0.5*i starts; "zero" uses zeros; "small_random"
#'   uses N(0, init_sd^2).
#' @param init_sd SD used when init_method = "small_random".
#' @param adapt Logical; adapt the coordinate-wise proposal SDs during warmup
#'   (Robbins-Monro, then frozen). Use a warmup of at least 1000 iterations.
#' @param target_accept Target acceptance rate for the adaptation.
#' @param penalize_intercept Logical; if FALSE (default) and the first column of Z
#'   is an intercept, it gets a vague N(0, intercept_sd^2) prior and is excluded
#'   from the Bayesian-lasso hierarchy. TRUE reproduces the previous behaviour.
#' @param intercept_sd Prior SD of the unpenalised intercept.
#' @param progress Logical; show progress for a single sequential chain.
#'
#' @note Optional entries of \code{hyperpar}: a_sigma and b_sigma (inverse-gamma
#'   prior of sigma^2; defaults 3 and 2).
#' @note \code{standardize = TRUE} centres and scales all columns of X and Z
#'   except the intercept (see scale.info() and scale.apply()); binary columns
#'   are included unless \code{scale.binary = FALSE}. \code{standardize = FALSE}
#'   leaves X and Z UNCHANGED. The MCMC runs on the scaled data (so the lasso
#'   penalty acts on the scaled coefficients), but fit.results, b.chains,
#'   beta.chains and lambda.chains are returned on the ORIGINAL covariate scale
#'   (back-transformed draw by draw). The incidence intercept then refers to
#'   z = 0 and the baseline hazards to x = 0 on the original scale. 
#'
#' @return An object of class "SMCM" with fit, posterior chains, posterior
#'   means, R-hat (for multiple chains), acceptance rates, model-fit results,
#'   processed data, priors and MCMC information.
#'
#' @example 
#' n <- 300
#' b.true <- c( 0.4, 0.5, 1)
#' beta.true <- c(1, 0.2 )
#' baseline.hazard.rates <- 1
#' intervals <- c(0, 16)
#' dat1 <- simSMCM(n, b.true, beta.true, baseline.hazard.rates, intervals, seed = 2026,
#'                 cens.start = 30/365, cens.end = 32,  prob.cov.X = 0.5, prob.cov.Z = 0.5, same.cov = TRUE )
#' 
#' priorPar <- list(r1 = 1, delta1 = 1, r2 = 1, delta2 = 1, a = 1, b = 1)
#' fit.smcm.args.J1 <- list(
#'   hyperpar = priorPar, nchains = 2L, nIter = 25000L, warmup = 5000L, thin = 20L,
#'   mcmc.parallel = "parLapply", standardize = TRUE, scale.binary = TRUE, save_loglik = 1L,
#'   init_method = "zero", penalize_intercept = FALSE,
#'   intercept_sd = 2,
#'   proposal_sd_b = 0.5, proposal_sd_beta = 0.5, proposal_sd_lambda = 0.5,
#'   probs = 0)
#' 
#' use_rcpp(TRUE)
#' system.time({
#'   fit <- do.call(fit.MCMC.SMCM, c(list(data = dat1), fit.smcm.args.J1,
#'                                   list(seed = 2026)) )
#' })
#' fit$fit.results
#' 
#' @references
#' Kızılaslan, F., Vitelli, V. (2025). Bayesian Semiparametric Mixture Cure (Frailty) Model,
#' arXiv \url{https://arxiv.org/abs/2512.08173}.
#'
#' @export
fit.MCMC.SMCM <- function(
    data,
    hyperpar = list(),
    nchains = 2L,
    nIter = 6000L,
    warmup = 1000L,
    thin = 20L,
    mcmc.parallel = "parLapply",
    standardize = TRUE,
    scale.binary =TRUE,
    probs,
    save_loglik = 0L,
    seed = NULL,
    store_tau = TRUE,
    proposal_sd_b = 1,
    proposal_sd_beta = 1,
    proposal_sd_lambda = 0.5,
    init_method = c("small_random", "zero", "legacy"),
    init_sd = 0.05,
    adapt = TRUE,
    target_accept = 0.35,
    penalize_intercept = FALSE,
    intercept_sd = 10,
    progress = interactive()
) {

  ## ------------------------------------------------------------------
  ## 0. Input checks
  ## ------------------------------------------------------------------
  init_method <- match.arg(init_method)
  mcmc.parallel <- match.arg(mcmc.parallel,
                             choices = c("lapply", "parLapply", "foreach"))

  nchains <- as.integer(nchains)
  nIter <- as.integer(nIter)
  warmup <- as.integer(warmup)
  thin <- as.integer(thin)
  save_loglik <- as.integer(isTRUE(as.logical(save_loglik)))

  if (nchains < 1L) stop("nchains must be >= 1.")
  if (nIter < 1L) stop("nIter must be >= 1.")
  if (warmup < 0L || warmup >= nIter) {
    stop("warmup must satisfy 0 <= warmup < nIter.")
  }
  if (thin < 1L) stop("thin must be >= 1.")
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
    stop("standardize must be TRUE or FALSE.")
  }

  required_data <- c("observed_time", "delta", "X", "Z")
  missing_data <- setdiff(required_data, names(data))
  if (length(missing_data)) {
    stop("Missing data component(s): ", paste(missing_data, collapse = ", "))
  }

  required_hyper <- c("r1", "delta1", "r2", "delta2", "a", "b")
  missing_hyper <- setdiff(required_hyper, names(hyperpar))
  if (length(missing_hyper)) {
    stop("Missing hyperparameter(s): ", paste(missing_hyper, collapse = ", "))
  }

  if (!exists("MCMC.SMCM.fast", mode = "function", inherits = TRUE)) {
    stop("MCMC.SMCM.fast() is not available. Source SMCM_fast.R first.")
  }

  if (!is.null(seed)) set.seed(seed)

  ## ------------------------------------------------------------------
  ## 1. Prepare survival data
  ## ------------------------------------------------------------------
  t <- as.numeric(data$observed_time)
  di <- as.integer(data$delta)
  X0 <- as.matrix(data$X)
  Z0 <- as.matrix(data$Z)

  n <- length(t)
  if (length(di) != n || nrow(X0) != n || nrow(Z0) != n) {
    stop("observed_time, delta, X and Z must have the same number of observations.")
  }
  if (any(!is.finite(t)) || any(t < 0)) {
    stop("observed_time must contain finite non-negative values.")
  }
  if (any(!di %in% c(0L, 1L))) {
    stop("delta must contain only 0 and 1.")
  }
  
  # standardize = TRUE (glmnet-style): centre/scale all columns except the intercept
  info.X <- info.Z <- NULL
  if (isTRUE(standardize)) {
    info.X <- scale.info(X0, scale.binary);  X <- scale.apply(X0, info.X)
    info.Z <- scale.info(Z0, scale.binary);  Z <- scale.apply(Z0, info.Z)
  } else {
    X <- X0
    Z <- Z0
  }
  
  survObj <- list(
    t = t,
    di = di,
    X = as.matrix(X),
    Z = as.matrix(Z)
  )

  p_beta <- ncol(survObj$X)
  p_b <- ncol(survObj$Z)

  ## ------------------------------------------------------------------
  ## 2. Piecewise-exponential cut points
  ## ------------------------------------------------------------------
  event.time.points <- sort(t[di == 1L])
  if (!length(event.time.points)) {
    stop("At least one observed event (delta = 1) is required to define cut points.")
  }

  probs <- as.numeric(probs)
  if (!length(probs) || any(!is.finite(probs))) {
    stop("probs must contain at least one finite value.")
  }
  if (any(probs < 0 | probs >= 1)) {
    stop("probs must lie in [0, 1).")
  }

  t.max <- max(t)
  lower.times <- t[t < t.max]
  if (length(lower.times)) {
    upper.bound <- 2 * t.max - max(lower.times)
  } else {
    upper.bound <- t.max + max(1, abs(t.max) * 0.01)
  }

  if (length(probs) == 1L && probs[1] == 0) {
    s <- c(0, upper.bound)
  } else {
    if (any(probs == 0)) {
      probs <- probs[probs > 0]
    }
    internal.cuts <- as.numeric(
      stats::quantile(event.time.points, probs = probs, names = FALSE)
    )
    internal.cuts <- internal.cuts[internal.cuts > 0 & internal.cuts < upper.bound]
    s <- sort(unique(c(0, internal.cuts, upper.bound)))
  }

  if (length(s) < 2L || any(diff(s) <= 0)) {
    stop("Could not construct strictly increasing piecewise-exponential cut points.")
  }

  JJ <- length(s) - 1L

  ## ------------------------------------------------------------------
  ## 3. Prior parameters
  ## ------------------------------------------------------------------
  priorPara <- list(
    s = s,
    J = JJ,
    r1 = hyperpar$r1,
    delta1 = hyperpar$delta1,
    r2 = hyperpar$r2,
    delta2 = hyperpar$delta2,
    a = hyperpar$a,
    b = hyperpar$b,
    a_sigma = if (is.null(hyperpar$a_sigma)) 3 else hyperpar$a_sigma,
    b_sigma = if (is.null(hyperpar$b_sigma)) 2 else hyperpar$b_sigma,
    penalize_intercept = penalize_intercept,
    intercept_sd = intercept_sd
  )

  ## ------------------------------------------------------------------
  ## 4. Base initial state
  ## ------------------------------------------------------------------
  initial <- list(
    beta.ini = rep(1, p_beta),
    b.ini = rep(1, p_b),
    sigmaSq = stats::runif(1L, 0.1, 10),
    sigmaStarSq = stats::runif(1L, 0.1, 10),
    etaSq = 1,
    etaStarSq = 1,
    lambda.ini = stats::rgamma(JJ, shape = 1, rate = 1)
  )

  initial$tauSq <- stats::rexp(p_b, rate = initial$etaSq / 2)
  initial$tauStarSq <- stats::rexp(p_beta, rate = initial$etaStarSq / 2)

  ## ------------------------------------------------------------------
  ## 5. Reproducible per-chain seeds
  ## ------------------------------------------------------------------
  if (is.null(seed)) {
    chain.seeds <- rep(NA_integer_, nchains)
  } else {
    chain.seeds <- sample.int(.Machine$integer.max, nchains)
  }

  ## ------------------------------------------------------------------
  ## 6. One chain runner -- used by every backend
  ## ------------------------------------------------------------------
  run_chain <- function(i, survObj, priorPara, initial,
                        nIter, warmup, thin, save_loglik, store_tau,
                        proposal_sd_b, proposal_sd_beta,
                        proposal_sd_lambda, init_method, init_sd,
                        chain.seeds, progress) {

    if (!is.na(chain.seeds[i])) set.seed(chain.seeds[i])

    initial.chain <- initial

    if (init_method == "legacy") {
      initial.chain$beta.ini <- rep(0.5 * i, ncol(survObj$X))
      initial.chain$b.ini <- rep(0.5 * i, ncol(survObj$Z))
    } else if (init_method == "zero") {
      initial.chain$beta.ini <- numeric(ncol(survObj$X))
      initial.chain$b.ini <- numeric(ncol(survObj$Z))
    } else {
      initial.chain$beta.ini <- stats::rnorm(ncol(survObj$X), 0, init_sd)
      initial.chain$b.ini <- stats::rnorm(ncol(survObj$Z), 0, init_sd)
    }

    ## chain-specific (overdispersed) starting values, drawn after the chain seed
    initial.chain$sigmaSq     <- stats::runif(1L, 0.1, 10)
    initial.chain$sigmaStarSq <- stats::runif(1L, 0.1, 10)
    initial.chain$lambda.ini  <- stats::rgamma(length(initial$lambda.ini), shape = 1, rate = 1)
    initial.chain$tauSq       <- stats::rexp(ncol(survObj$Z), rate = initial.chain$etaSq / 2)
    initial.chain$tauStarSq   <- stats::rexp(ncol(survObj$X), rate = initial.chain$etaStarSq / 2)

    MCMC.SMCM.fast(
      survObj = survObj,
      priorPara = priorPara,
      initial = initial.chain,
      nIter = nIter,
      warmup = warmup,
      thin = thin,
      save_loglik = save_loglik,
      store_tau = store_tau,
      proposal_sd_b = proposal_sd_b,
      proposal_sd_beta = proposal_sd_beta,
      proposal_sd_lambda = proposal_sd_lambda,
      adapt = adapt,                    # taken from the enclosing function environment
      target_accept = target_accept,
      progress = progress
    )
  }

  ## ------------------------------------------------------------------
  ## Helper: export only the functions needed by MCMC.SMCM.fast to PSOCK
  ## workers. This is used by both parLapply and foreach.
  ## ------------------------------------------------------------------
  setup_fast_workers <- function(cl) {
    fast.function.names <- c(
      "MCMC.SMCM.fast",
      "cumulative_hazard_matrix", "smcm_hazard_from_state",
      "softplus", "softplus_smcm",
      "loglik_b_part", "loglik_beta_part", "loglik_lambda_part",
      "smcm_loglik_b", "smcm_loglik_beta", "smcm_loglik_lambda",
      "logLik_SMCM_fast",
      "UpdateRP.b.fast", "UpdateRP.beta.fast", "UpdateBH.lambda.fast",
      "UpdateTau.fast", "UpdateTauStar.fast",
      "UpdateSigma.fast", "UpdateSigmaStar.fast",
      "UpdateEta.fast", "UpdateEtaStar.fast"
    )

    search.env <- environment()
    available <- fast.function.names[
      vapply(fast.function.names,
             function(z) exists(z, envir = search.env,
                                mode = "function", inherits = TRUE),
             logical(1))
    ]

    worker.env <- new.env(parent = emptyenv())
    for (nm in available) {
      assign(nm,
             get(nm, envir = search.env, mode = "function", inherits = TRUE),
             envir = worker.env)
    }
    if (length(available)) {
      parallel::clusterExport(cl, available, envir = worker.env)
      
      if (isTRUE(getOption("smcm.use.rcpp"))) {
        compile_on_worker <- function(i, f) {
          cd <- file.path(dirname(f), ".smcm_cache", i)       # one cache folder per worker, reused across fits
          dir.create(cd, recursive = TRUE, showWarnings = FALSE)
          Rcpp::sourceCpp(f, env = .GlobalEnv, cacheDir = cd)
          invisible(TRUE)
        }
        environment(compile_on_worker) <- globalenv()         # do not ship the enclosing environment
        parallel::clusterApply(cl, seq_along(cl), compile_on_worker, getOption("smcm.cpp.file"))
      }
    }
    

    invisible(available)
  }

  ## ------------------------------------------------------------------
  ## 7. Run chains
  ## ------------------------------------------------------------------
  start.mcmc.smc.time <- Sys.time()

  if (nchains == 1L) {
    fit.list <- list(
      run_chain(
        1L, survObj, priorPara, initial,
        nIter, warmup, thin, save_loglik, store_tau,
        proposal_sd_b, proposal_sd_beta, proposal_sd_lambda,
        init_method, init_sd, chain.seeds,
        progress = isTRUE(progress)
      )
    )

  } else if (mcmc.parallel == "lapply") {
    fit.list <- lapply(
      seq_len(nchains),
      run_chain,
      survObj = survObj,
      priorPara = priorPara,
      initial = initial,
      nIter = nIter,
      warmup = warmup,
      thin = thin,
      save_loglik = save_loglik,
      store_tau = store_tau,
      proposal_sd_b = proposal_sd_b,
      proposal_sd_beta = proposal_sd_beta,
      proposal_sd_lambda = proposal_sd_lambda,
      init_method = init_method,
      init_sd = init_sd,
      chain.seeds = chain.seeds,
      progress = FALSE
    )

  } else if (mcmc.parallel == "parLapply") {
    cl <- parallel::makeCluster(nchains)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    setup_fast_workers(cl)

    fit.list <- parallel::parLapply(
      cl,
      seq_len(nchains),
      run_chain,
      survObj = survObj,
      priorPara = priorPara,
      initial = initial,
      nIter = nIter,
      warmup = warmup,
      thin = thin,
      save_loglik = save_loglik,
      store_tau = store_tau,
      proposal_sd_b = proposal_sd_b,
      proposal_sd_beta = proposal_sd_beta,
      proposal_sd_lambda = proposal_sd_lambda,
      init_method = init_method,
      init_sd = init_sd,
      chain.seeds = chain.seeds,
      progress = FALSE
    )

  } else { # foreach
    if (!requireNamespace("foreach", quietly = TRUE) ||
        !requireNamespace("doParallel", quietly = TRUE)) {
      stop("Packages 'foreach' and 'doParallel' are required for mcmc.parallel='foreach'.")
    }

    cl <- parallel::makeCluster(nchains)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    setup_fast_workers(cl)
    doParallel::registerDoParallel(cl)

    `%dopar%` <- foreach::`%dopar%`
    fit.list <- foreach::foreach(i = seq_len(nchains)) %dopar% {
      run_chain(
        i, survObj, priorPara, initial,
        nIter, warmup, thin, save_loglik, store_tau,
        proposal_sd_b, proposal_sd_beta, proposal_sd_lambda,
        init_method, init_sd, chain.seeds,
        progress = FALSE
      )
    }
  }

  mcmc.smc.time <- as.numeric(
    difftime(Sys.time(), start.mcmc.smc.time, units = "mins")
  )

  ## ------------------------------------------------------------------
  ## 8. Extract saved chains
  ## ------------------------------------------------------------------
  b.list <- lapply(fit.list, function(x) as.matrix(x$b.p))
  beta.list <- lapply(fit.list, function(x) as.matrix(x$beta.p))
  lambda.list <- lapply(fit.list, function(x) as.matrix(x$lambda.p))

  if (!length(b.list[[1]]) || nrow(b.list[[1]]) < 1L) {
    stop("No posterior draws were saved. Check nIter, warmup and thin.")
  }

  b.hat.mcmc <- colMeans(do.call(rbind, b.list))
  beta.hat.mcmc <- colMeans(do.call(rbind, beta.list))
  lambda.hat.mcmc <- colMeans(do.call(rbind, lambda.list))

  ## ------------------------------------------------------------------
  ## 9. R-hat from already processed chains
  ## ------------------------------------------------------------------
  rhat_one_block <- function(chain.list, prefix) {
    p <- ncol(chain.list[[1]])

    if (length(chain.list) < 2L) {
      out <- matrix(NA_real_, nrow = 1L, ncol = p)
      colnames(out) <- paste0(prefix, "[", seq_len(p), "]")
      return(out)
    }

    n.draw <- vapply(chain.list, nrow, integer(1))
    if (length(unique(n.draw)) != 1L || n.draw[1] < 2L) {
      out <- matrix(NA_real_, nrow = 1L, ncol = p)
      colnames(out) <- paste0(prefix, "[", seq_len(p), "]")
      return(out)
    }

    if (requireNamespace("coda", quietly = TRUE)) {
      ans <- tryCatch({
        ml <- coda::mcmc.list(
          lapply(chain.list, function(z) coda::mcmc(as.matrix(z)))
        )
        psrf <- coda::gelman.diag(
          ml, autoburnin = FALSE, multivariate = FALSE
        )$psrf[, 1]
        as.numeric(psrf)
      }, error = function(e) rep(NA_real_, p))
    } else {
      # Classical Gelman-Rubin calculation as a dependency-free fallback.
      m <- length(chain.list)
      nn <- n.draw[1]
      chain.means <- do.call(rbind, lapply(chain.list, colMeans))
      chain.vars <- do.call(rbind, lapply(chain.list, function(z) {
        apply(z, 2L, stats::var)
      }))
      W <- colMeans(chain.vars)
      B <- nn * apply(chain.means, 2L, stats::var)
      var.plus <- ((nn - 1) / nn) * W + B / nn
      ans <- sqrt(var.plus / W)
      ans[!is.finite(ans)] <- NA_real_
    }

    out <- matrix(ans, nrow = 1L)
    colnames(out) <- paste0(prefix, "[", seq_len(p), "]")
    out
  }

  if (nchains > 1L) {
    hyper.list <- lapply(fit.list, function(x) cbind(
      sigmaSq = x$mcmcOutcome$sigmaSq.p, sigmaStarSq = x$mcmcOutcome$sigmaStarSq.p,
      etaSq = x$mcmcOutcome$etaSq.p, etaStarSq = x$mcmcOutcome$etaStarSq.p))
    rhat.hyper <- rhat_one_block(hyper.list, "hyper")
    colnames(rhat.hyper) <- colnames(hyper.list[[1]])

    Rhat <- cbind(
      rhat_one_block(b.list, "b"),
      rhat_one_block(beta.list, "beta"),
      rhat_one_block(lambda.list, "lambda"),
      rhat.hyper
    )
  } else {
    Rhat <- NULL
  }

  ## Effective sample sizes (summed over chains); requires coda
  ess <- NULL
  if (requireNamespace("coda", quietly = TRUE)) {
    ess_fun <- function(chain.list) {
      tryCatch(
        coda::effectiveSize(coda::mcmc.list(lapply(chain.list, function(z) coda::mcmc(as.matrix(z))))),
        error = function(e) NULL)
    }
    ess <- list(b = ess_fun(b.list), beta = ess_fun(beta.list), lambda = ess_fun(lambda.list))
  }

  ## ------------------------------------------------------------------
  ## 10. Observation-level log-likelihood and DIC/LPML
  ## ------------------------------------------------------------------
  loo.results <- NULL
  loglik.list <- NULL

  if (save_loglik == 1L) {
    loglik.list <- lapply(fit.list, function(x) as.matrix(x$loglik.matrix))

    if (!exists("mcmc.dic.lpml", mode = "function", inherits = TRUE)) {
      warning("mcmc.dic.lpml() is not available; loo.results set to NULL.")
    } else {
      loglik.for.summary <- if (nchains == 1L) loglik.list[[1]] else loglik.list
      loo.results <- mcmc.dic.lpml(
        survObj,
        priorPara,
        loglik.for.summary,
        b.hat.mcmc,
        beta.hat.mcmc,
        lambda.hat.mcmc,
        stan.model = FALSE,
        frailty = FALSE
      )
    }
  }

  ## ------------------------------------------------------------------
  ## 11. Acceptance-rate summaries
  ## ------------------------------------------------------------------
  acceptance.b <- do.call(
    rbind,
    lapply(fit.list, function(x) as.numeric(x$mcmcOutcome$acceptance.b))
  )
  acceptance.beta <- do.call(
    rbind,
    lapply(fit.list, function(x) as.numeric(x$mcmcOutcome$acceptance.beta))
  )
  acceptance.lambda <- do.call(
    rbind,
    lapply(fit.list, function(x) as.numeric(x$mcmcOutcome$acceptance.lambda))
  )

  rownames(acceptance.b) <- paste0("chain", seq_len(nchains))
  rownames(acceptance.beta) <- paste0("chain", seq_len(nchains))
  rownames(acceptance.lambda) <- paste0("chain", seq_len(nchains))

  colnames(acceptance.b) <- paste0("b[", seq_len(ncol(acceptance.b)), "]")
  colnames(acceptance.beta) <- paste0("beta[", seq_len(ncol(acceptance.beta)), "]")
  colnames(acceptance.lambda) <- paste0("lambda[", seq_len(ncol(acceptance.lambda)), "]")

  ## ------------------------------------------------------------------
  ## 12. Posterior incidence probability and final result
  ## ------------------------------------------------------------------
  # Posterior mean of pi(z_i) = plogis(z_i' b), averaged over all saved draws
  b.draws <- do.call(rbind, b.list)
  pz.hat.mcmc <- rowMeans(stats::plogis(survObj$Z %*% t(b.draws)))

  ## Back-transform to the ORIGINAL covariate scale (only if the fit used scaled X, Z).
  scaling <- list(applied = isTRUE(standardize), X = info.X, Z = info.Z)
  if (scaling$applied) {
    cs <- scaling$X$center / scaling$X$scale
    # baseline hazard absorbs the centring constant: lambda_orig = lambda_s * exp(-sum(m_j beta_j,s / s_j))
    lambda.out <- Map(function(L, Bt) L * exp(-drop(Bt %*% cs)), lambda.list, beta.list)
    b.out      <- lapply(b.list,    backtransform.b,    info = scaling$Z)
    beta.out   <- lapply(beta.list, backtransform.beta, info = scaling$X)
  } else {
    b.out <- b.list;  beta.out <- beta.list;  lambda.out <- lambda.list
  }
  b.hat.out      <- colMeans(do.call(rbind, b.out))
  beta.hat.out   <- colMeans(do.call(rbind, beta.out))
  lambda.hat.out <- colMeans(do.call(rbind, lambda.out))

  fit.results <- c(
    b.hat = b.hat.out,
    beta.hat = beta.hat.out,
    lambda.hat = lambda.hat.out,
    pz.hat = mean(pz.hat.mcmc),
    time = mcmc.smc.time
  )

  one <- function(lst) if (nchains == 1L) lst[[1]] else lst # as in the fit_MCMC_SMCFM.R1
  ret <- list(
    fit = one(fit.list), 
    fit.results = fit.results,
    b.chains = one(b.out), 
    beta.chains = one(beta.out), 
    lambda.chains = one(lambda.out), 
    Rhat = Rhat,
    ess = ess,
    scaling = scaling,
    loo.results = loo.results,
    loglik.matrix.chains =  if (save_loglik == 1L) one(loglik.list) else NULL, #loglik.out,
    acceptance = list(
      b = acceptance.b,
      beta = acceptance.beta,
      lambda = acceptance.lambda
    ),
    data = data,#survObj,
    priors = priorPara,
    mcmc.info = list(
      nchains = nchains,
      nIter = nIter,
      warmup = warmup,
      thin = thin,
      n_saved_per_chain = nrow(b.list[[1]]),
      parallel = mcmc.parallel,
      runtime_minutes = mcmc.smc.time,
      seed = seed,
      chain_seeds = chain.seeds,
      init_method = init_method,
      proposal_sd_b = proposal_sd_b,
      proposal_sd_beta = proposal_sd_beta,
      proposal_sd_lambda = proposal_sd_lambda,
      store_tau = store_tau,
      save_loglik = save_loglik,
      adapt = adapt,
      target_accept = target_accept,
      penalize_intercept = penalize_intercept,
      intercept_sd = intercept_sd,
      proposal_sd_final = lapply(fit.list, function(x) x$mcmcOutcome$proposal.sd)
    )
  )

  class(ret) <- "SMCM"
  ret
}
## Helpers to return posterior draws to the ORIGINAL covariate scale when the model
## was fitted on scale.apply scaled X and Z (standardize = TRUE).
## Columns with non-integer values are centred/scaled; integer columns (binary
## dummies, intercept) are untouched, i.e. centre 0 and scale 1.
#' @noRd
scale.info <- function(M0, scale.binary = TRUE) {
  sdv   <- apply(M0, 2, stats::sd)
  const <- sdv == 0                                          # intercept / constant columns: never scaled
  skip  <- if (scale.binary) const else const | apply(M0, 2, function(x) all(x == floor(x)))
  list(center = ifelse(skip, 0, colMeans(M0)), scale = ifelse(skip, 1, sdv))
}
#' @noRd
scale.apply <- function(M0, info) sweep(sweep(M0, 2, info$center, "-"), 2, info$scale, "/")
# scale.info <- function(M0) {
#   cont <- apply(M0, 2, function(x) !all(x == floor(x)))
#   list(center = ifelse(cont, colMeans(M0), 0),
#        scale  = ifelse(cont, apply(M0, 2, stats::sd), 1))
# }
#' @noRd
backtransform.b <- function(B, info) {          # B: draws x p_b, first column = intercept
  s <- info$scale; m <- info$center
  out <- sweep(B, 2, s, "/")                    # slopes: b_j / s_j
  out[, 1] <- B[, 1] - as.vector(B[, -1, drop = FALSE] %*% (m[-1] / s[-1]))  # b0 - sum(m_j b_j / s_j)
  out
}
#' @noRd
backtransform.beta <- function(Bt, info) sweep(Bt, 2, info$scale, "/")       # beta_j / s_j
