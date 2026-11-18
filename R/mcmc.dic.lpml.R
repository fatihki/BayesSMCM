#' mcmc.dic.lpml
#' @title DIC, LPML, LOO-CV, AIC and BIC from MCMC output (SMCM and SMCFM)
#'
#' @param survObj list with t, di, X, Z (the data the MCMC was run on, i.e. SCALED if standardize = TRUE).
#' @param priorPara list with the piecewise-exponential cut points s (as used in the fit).
#' @param loglik.matrix.chains observation-level log-likelihood draws: EITHER a list with one
#'   (draws x n) matrix per chain, OR a single (draws x n) matrix (one chain, or the already combined
#'   chains of a Stan fit).
#' @param b.hat.mcmc,beta.hat.mcmc,lambda.hat.mcmc posterior means (same scale as survObj).
#' @param stan.model kept for backward compatibility; a list is stacked, a matrix is used as is.
#' @param frailty FALSE -> SMCM (logLik_SMCM), TRUE -> SMCFM (logLik_SMCFM, needs \code{theta.hat.mcmc}).
#' @param ... \code{theta.hat.mcmc} (required when frailty = TRUE); optionally \code{chain_id}
#'   (integer vector, one entry per row of a combined matrix) used for the relative ESS in loo.
#'   
#' @return list(loo_result, lpml, cpo, LL, LL_pm, k, AIC, BIC, AIC_avg, BIC_avg, p_D, DIC,
#'   p_V, DIC_V, Dev_average, Dev_at_pm). Smaller DIC / AIC / BIC / looic and larger LPML are better.
#' @noRd
mcmc.dic.lpml <- function(survObj, priorPara, loglik.matrix.chains,
                          b.hat.mcmc, beta.hat.mcmc, lambda.hat.mcmc,
                          stan.model = FALSE, frailty = FALSE, ...) {
  args <- list(...)

  theta.hat.mcmc <- NA_real_
  if (isTRUE(frailty)) {
    theta.hat.mcmc <- args$theta.hat.mcmc
    if (is.null(theta.hat.mcmc) || length(theta.hat.mcmc) != 1L || !is.finite(theta.hat.mcmc)) {
      stop("frailty = TRUE requires a finite scalar 'theta.hat.mcmc'.")
    }
  }

  ## ---- stack the chains ------------------------------------------------------------------
  if (is.list(loglik.matrix.chains) && !is.data.frame(loglik.matrix.chains)) {
    mats          <- lapply(loglik.matrix.chains, as.matrix)
    loglik.matrix <- do.call(rbind, mats)
    chain_id      <- rep(seq_along(mats), vapply(mats, nrow, integer(1)))
  } else {
    loglik.matrix <- as.matrix(loglik.matrix.chains)
    chain_id      <- if (!is.null(args$chain_id)) as.integer(args$chain_id) else rep(1L, nrow(loglik.matrix))
  }
  if (length(chain_id) != nrow(loglik.matrix)) stop("chain_id must have one entry per row of the log-likelihood matrix.")
  if (ncol(loglik.matrix) != length(survObj$t)) {
    stop("The log-likelihood matrix must have one column per observation (", length(survObj$t), ").")
  }
  if (any(!is.finite(loglik.matrix))) stop("Non-finite values in the log-likelihood matrix.")
  S <- nrow(loglik.matrix)

  ## ---- CPO and LPML (harmonic-mean estimator, on the log scale) --------------------------
  # CPO_i = 1 / mean_s( 1 / L_is )  =>  log CPO_i = -( logsumexp_s(-ll_is) - log S )
  logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
  log_cpo   <- apply(loglik.matrix, 2L, function(x) -(logsumexp(-x) - log(length(x))))
  lpml      <- sum(log_cpo)

  ## ---- PSIS-LOO with the relative ESS ----------------------------------------------------
  r_eff <- tryCatch({
    L <- exp(sweep(loglik.matrix, 2L, apply(loglik.matrix, 2L, max), "-"))   # scaling by a constant per column does not change r_eff
    loo::relative_eff(L, chain_id = chain_id, cores = 1)
  }, error = function(e) NA)
  loo_result <- loo::loo(loglik.matrix, r_eff = r_eff, cores = 1)

  ## ---- deviance based criteria -----------------------------------------------------------
  Dev_draws   <- -2 * rowSums(loglik.matrix)
  Dev_average <- mean(Dev_draws)                                     # posterior mean deviance
  LL          <- mean(rowSums(loglik.matrix))                        # posterior mean log-likelihood

  Dev_at_pm <- if (!isTRUE(frailty)) {
    -2 * sum(logLik_SMCM_fast(survObj, priorPara, b.coef = b.hat.mcmc, beta.coef = beta.hat.mcmc,
                         lambda.coef = lambda.hat.mcmc))
  } else {
    -2 * sum(logLik_SMCFM(survObj, priorPara, b.coef = b.hat.mcmc, beta.coef = beta.hat.mcmc,
                          lambda.coef = lambda.hat.mcmc, theta.coef = theta.hat.mcmc))
  }
  LL_pm <- -Dev_at_pm / 2                                            # log-likelihood at the posterior means

  p_D <- Dev_average - Dev_at_pm
  DIC <- Dev_at_pm + 2 * p_D                                         # = Dev_average + p_D
  p_V   <- stats::var(Dev_draws) / 2
  DIC_V <- Dev_average + p_V

  ## ---- AIC / BIC (plug-in at the posterior means) ----------------------------------------
  k <- length(b.hat.mcmc) + length(beta.hat.mcmc) + length(lambda.hat.mcmc) + ifelse(isTRUE(frailty), 1L, 0L)
  n <- length(survObj$t)
  AIC <- 2 * k - 2 * LL_pm
  BIC <- k * log(n) - 2 * LL_pm


  list(loo_result = loo_result, lpml = lpml, cpo = exp(log_cpo),
       LL = LL, LL_pm = LL_pm, k = k,
       AIC = AIC, BIC = BIC,
       p_D = p_D, DIC = DIC, p_V = p_V, DIC_V = DIC_V,
       Dev_average = Dev_average, Dev_at_pm = Dev_at_pm)
}
