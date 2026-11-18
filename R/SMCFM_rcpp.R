## Rcpp switch for the SMCFM (mixture cure FRAILTY) sampler.
## Source this AFTER SMCM_fast_R1.R and SMCFM_fast_R1.R (those files define the UpdateRP.* /
## UpdateBH.* / UpdateTheta.* functions and would overwrite the Rcpp versions again).
##
##   source("SMCM_fast_R1.R")
##   source("SMCFM_fast_R1.R")
##   source("SMCM_rcpp.R")      # optional: Rcpp versions of the SMCM updates AND of UpdateTau.fast /
##                              #           UpdateTauStar.fast, which the SMCFM sampler also uses
##   source("SMCFM_rcpp.R")     # compiles SMCFM_rcpp.cpp and calls use_rcpp_smcfm(TRUE)
##   use_rcpp_smcfm(FALSE)      # back to the pure-R SMCFM updates (e.g. for comparisons)
##
## The lasso hierarchy (UpdateTau*, UpdateSigma*, UpdateEta*) is shared with SMCM and is not
## touched here; whatever SMCM_rcpp.R / use_rcpp() set for those functions is used by both models.
##
## Needs: install.packages("Rcpp") and a C++ compiler (Rtools on Windows,
## Xcode command line tools on macOS, build-essential on Linux).

stopifnot(exists("UpdateRP.b.SMCFM.fast"), exists("UpdateRP.beta.SMCFM.fast"),
          exists("UpdateBH.lambda.SMCFM.fast"), exists("UpdateTheta.SMCFM.fast"))
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Please install the Rcpp package.")

smcfm.cpp.file <- normalizePath("SMCFM_rcpp.cpp", mustWork = TRUE)
Rcpp::sourceCpp(smcfm.cpp.file)           # first call compiles (about 30 s); later calls reuse the cache

## keep the original R versions.  The Rcpp wrappers below carry the attribute "smcfm.rcpp", so that
## sourcing this file AGAIN (without re-sourcing SMCFM_fast_R1.R) never stores a wrapper as the "pure R"
## version: only genuine R functions are (re)captured, otherwise the copy kept from the first sourcing is used.
.smcfm.names <- c("UpdateRP.b.SMCFM.fast", "UpdateRP.beta.SMCFM.fast",
                  "UpdateBH.lambda.SMCFM.fast", "UpdateTheta.SMCFM.fast")
if (!exists(".smcfm_r", inherits = FALSE)) .smcfm_r <- list()
for (nm in .smcfm.names) {
  cur <- get(nm, envir = .GlobalEnv, mode = "function")
  if (!isTRUE(attr(cur, "smcfm.rcpp"))) {
    .smcfm_r[[nm]] <- cur                                    # genuine R version: (re)capture it
  } else if (is.null(.smcfm_r[[nm]])) {
    stop(nm, " is already an Rcpp wrapper and no R original is stored: ",
         "source SMCFM_fast_R1.R again before sourcing SMCFM_rcpp.R.")
  }                                                          # else keep the stored original
}

## drop-in replacements with the same arguments and return values as the R versions
.smcfm_cpp <- list(
  UpdateRP.b.SMCFM.fast = function(Z, b, Zb, logS, event, sd.b, proposal.sd = rep(1, length(b))) {
    if (length(proposal.sd) == 1L) proposal.sd <- rep(proposal.sd, length(b))
    if (!is.double(Z)) storage.mode(Z) <- "double"            # avoids a copy on every call when already double
    UpdateRP_b_SMCFM_cpp(Z, b, Zb, logS, as.integer(event), sd.b, proposal.sd)
  },
  UpdateRP.beta.SMCFM.fast = function(X, beta, Xbeta, Zb, H0, theta, event, sd.beta,
                                      proposal.sd = rep(1, length(beta))) {
    if (length(proposal.sd) == 1L) proposal.sd <- rep(proposal.sd, length(beta))
    if (!is.double(X)) storage.mode(X) <- "double"
    UpdateRP_beta_SMCFM_cpp(X, beta, Xbeta, Zb, H0, as.numeric(theta), as.integer(event),
                            sd.beta, proposal.sd)
  },
  UpdateBH.lambda.SMCFM.fast = function(lambda, H0, h0, hazard_info, Xbeta, Zb, theta, event,
                                        a, b, proposal.sd = rep(1, length(lambda))) {
    J <- length(lambda)
    if (length(proposal.sd) == 1L) proposal.sd <- rep(proposal.sd, J)
    if (length(proposal.sd) != J) stop("proposal.sd must have length 1 or length(lambda).")
    if (any(!is.finite(proposal.sd)) || any(proposal.sd <= 0)) stop("proposal.sd must contain positive finite values.")
    if (any(!is.finite(lambda)) || any(lambda <= 0)) stop("lambda must contain positive finite values.")
    UpdateBH_lambda_SMCFM_cpp(lambda, H0, h0, hazard_info$A, as.integer(hazard_info$interval_id),
                              Xbeta, Zb, as.numeric(theta), as.integer(event), a, b, proposal.sd)
  },
  UpdateTheta.SMCFM.fast = function(theta, H1, Zb, event, c.shape, d.rate, proposal.sd = 0.5)
    UpdateTheta_SMCFM_cpp(as.numeric(theta), H1, Zb, as.integer(event), c.shape, d.rate, proposal.sd)
)

.smcfm_cpp <- lapply(.smcfm_cpp, function(f) { attr(f, "smcfm.rcpp") <- TRUE; f })   # tag the wrappers

use_rcpp_smcfm <- function(on = TRUE) {
  fns <- if (on) .smcfm_cpp else .smcfm_r
  for (nm in names(fns)) assign(nm, fns[[nm]], envir = .GlobalEnv)
  options(smcfm.use.rcpp = on, smcfm.cpp.file = smcfm.cpp.file)
  message("SMCFM updates: ", if (on) "Rcpp" else "pure R")
  invisible(on)
}

use_rcpp_smcfm(TRUE)


## ---------------------------------------------------------------------------
## Parallel chains (mcmc.parallel = "parLapply" or "foreach"):
## compiled functions cannot be sent to PSOCK workers, so each worker must compile
## the C++ file(s) itself.  In fit_MCMC_SMCFM_R1.R, inside setup_fast_workers(), add this
## right after the parallel::clusterExport(cl, available, envir = worker.env) call:
##
##   cpp.files <- c(if (isTRUE(getOption("smcm.use.rcpp")))   getOption("smcm.cpp.file"),    # UpdateTau_cpp (SMCM_rcpp.cpp)
##                  if (isTRUE(getOption("smcfm.use.rcpp")))  getOption("smcfm.cpp.file"))   # SMCFM updates
##   if (length(cpp.files)) {
##     compile_on_worker <- function(i, files) {
##       for (f in files) {
##         cd <- file.path(dirname(f), ".smcm_cache", tools::file_path_sans_ext(basename(f)), i)
##         dir.create(cd, recursive = TRUE, showWarnings = FALSE)   # one cache folder per worker and file,
##         Rcpp::sourceCpp(f, env = .GlobalEnv, cacheDir = cd)      # reused across fits AND R sessions
##       }
##       invisible(TRUE)
##     }
##     environment(compile_on_worker) <- globalenv()                # do not ship the enclosing environment
##     parallel::clusterApply(cl, seq_along(cl), compile_on_worker, cpp.files)
##   }
##
## (the patched file is provided as fit_MCMC_SMCFM_R1.R).  Until you use it, run with
## mcmc.parallel = "lapply" (chains run one after another).
## ---------------------------------------------------------------------------