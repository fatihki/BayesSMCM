# # R/zzz.R
# .onLoad <- function(libname, pkgname) {
#   # Enable automatic caching of compiled Stan models
#   rstan::rstan_options(auto_write = TRUE)
# 
#   # Enable parallel sampling using all available cores
#   options(mc.cores = parallel::detectCores())
# }


# R/zzz.R

.onLoad <- function(libname, pkgname) {
  # -----------------------------
  # 1. rstan setup
  # -----------------------------
  # # Enable automatic caching of compiled Stan models
  # rstan::rstan_options(auto_write = TRUE)
  # 
  # # Enable parallel sampling using all available cores
  # options(mc.cores = parallel::detectCores())
  
  # -----------------------------
  # 2. Suppress R CMD check warnings for base functions
  # -----------------------------
  if (getRversion() >= "2.15.1") {
    utils::globalVariables(
      c( 
        "runif", "rexp", "rgamma", "dgamma", "dnorm","rnorm", "rweibull", "pweibull", "rbinom", 
        "quantile", "scale", "mean", "sd", "window"
      )
    )
  }
}
