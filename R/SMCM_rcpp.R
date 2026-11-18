## Rcpp switch for the SMCM sampler.
## Source this AFTER SMCM_fast_R1.R (and after any other file that defines the
## UpdateRP.* / UpdateBH.* / UpdateTau* functions), because those files would
## overwrite the Rcpp versions again.
##
##   source("SMCM_fast_R1.R")
##   source("SMCM_rcpp.R")      # compiles SMCM_rcpp.cpp and calls use_rcpp(TRUE)
##   use_rcpp(FALSE)            # back to the pure-R updates (e.g. for comparisons)
##
## Needs: install.packages("Rcpp") and a C++ compiler (Rtools on Windows,
## Xcode command line tools on macOS, build-essential on Linux).

stopifnot(exists("UpdateRP.b.fast"), exists("UpdateRP.beta.fast"), exists("UpdateBH.lambda.fast"),
          exists("UpdateTau.fast"), exists("UpdateTauStar.fast"))
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Please install the Rcpp package.")

smcm.cpp.file <- normalizePath("SMCM_rcpp.cpp", mustWork = TRUE)
Rcpp::sourceCpp(smcm.cpp.file)            # first call compiles (about 30 s); later calls reuse the cache

## keep the original R versions
.smcm_r <- list(
  UpdateRP.b.fast      = UpdateRP.b.fast,
  UpdateRP.beta.fast   = UpdateRP.beta.fast,
  UpdateBH.lambda.fast = UpdateBH.lambda.fast,
  UpdateTau.fast       = UpdateTau.fast,
  UpdateTauStar.fast   = UpdateTauStar.fast
)

## drop-in replacements with the same arguments and return values as the R versions
.smcm_cpp <- list(
  UpdateRP.b.fast = function(Z, b, Zb, H1, delta, sd.b, proposal.sd = rep(1, length(b))) {
    if (length(proposal.sd) == 1L) proposal.sd <- rep(proposal.sd, length(b))
    if (!is.double(Z)) storage.mode(Z) <- "double"            # avoids a copy on every call when already double
    UpdateRP_b_cpp(Z, b, Zb, H1, as.integer(delta), sd.b, proposal.sd)
  },
  UpdateRP.beta.fast = function(X, beta, Xbeta, Zb, H0, delta, sd.beta,
                                proposal.sd = rep(1, length(beta))) {
    if (length(proposal.sd) == 1L) proposal.sd <- rep(proposal.sd, length(beta))
    if (!is.double(X)) storage.mode(X) <- "double"
    UpdateRP_beta_cpp(X, beta, Xbeta, Zb, H0, as.integer(delta), sd.beta, proposal.sd)
  },
  UpdateBH.lambda.fast = function(lambda, H0, h0, hazard_info, Xbeta, Zb, delta,
                                  a, b, proposal.sd = rep(1, length(lambda))) {
    if (length(proposal.sd) == 1L) proposal.sd <- rep(proposal.sd, length(lambda))
    if (any(!is.finite(lambda)) || any(lambda <= 0)) stop("lambda must contain positive finite values.")
    UpdateBH_lambda_cpp(lambda, H0, h0, hazard_info$A, as.integer(hazard_info$interval_id),
                        Xbeta, Zb, as.integer(delta), a, b, proposal.sd)
  },
  UpdateTau.fast = function(b, sigmaSq, etaSq, eps = 1e-12)
    UpdateTau_cpp(b, sigmaSq, etaSq, eps),
  UpdateTauStar.fast = function(beta, sigmaStarSq, etaStarSq, eps = 1e-12)
    UpdateTau_cpp(beta, sigmaStarSq, etaStarSq, eps)
)

use_rcpp <- function(on = TRUE) {
  fns <- if (on) .smcm_cpp else .smcm_r
  for (nm in names(fns)) assign(nm, fns[[nm]], envir = .GlobalEnv)
  options(smcm.use.rcpp = on, smcm.cpp.file = smcm.cpp.file)   
  message("SMCM updates: ", if (on) "Rcpp" else "pure R")
  invisible(on)
}

use_rcpp(TRUE)


## ---------------------------------------------------------------------------
## Parallel chains (mcmc.parallel = "parLapply" or "foreach"):
## compiled functions cannot be sent to PSOCK workers, so each worker must compile
## the C++ file itself. Add this inside setup_fast_workers() in fit_MCMC_SMCM_R1.R,
## right after the parallel::clusterExport(...) call:
##
##   if (isTRUE(getOption("smcm.use.rcpp"))) {
##     parallel::clusterCall(cl, function(f) {
##       Rcpp::sourceCpp(f, env = .GlobalEnv, cacheDir = tempfile("smcm_"))
##       invisible(TRUE)
##     }, getOption("smcm.cpp.file"))
##   }
##
## Until you add it, use mcmc.parallel = "lapply" (chains run one after another).
## ---------------------------------------------------------------------------
