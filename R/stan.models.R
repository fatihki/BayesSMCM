#' @title Internal Pre-compiled Stan models for SMCM and SMCFM
#' @description Compiles the Stan model once and reuses it in package functions.
#' @keywords internal
#'
#' @noRd
SMCM.model            <- rstan::stan_model(
  file = system.file("stan", "SMCM.stan", package = "BayesSMCM"))
#' @noRd
SMCFM.model            <- rstan::stan_model(
  file = system.file("stan", "SMCFM.stan", package = "BayesSMCM"))
#' @noRd
HSMCM.model           <- rstan::stan_model(
  file = system.file("stan", "HSMCM.stan", package = "BayesSMCM") ) 
#' @noRd
HSMCFM.model <- rstan::stan_model(
  file = system.file("stan", "HSMCFM.stan", package = "BayesSMCM")
)
