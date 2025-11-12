#' @title Highest Posterior Density (HPD) intervals for all MCMC chains of SMCM and SMCFM.
#' 
#' @description
#' It is calculated by using \pkg{coda} via \link[coda]{HPDinterval}.
#' 
#' @name HPD.interval.mcmc
#' 
#' @importFrom coda mcmc as.mcmc mcmc.list HPDinterval
#' 
#' @param fit.mcmc.chains MCMC results for all the chains based on our MCMC results.
#' @param pars Parameters can be "beta.p" or "b.p" or "lambda.p".
#' @param warmup Number of warmup iterations.
#' @param thin Thinning interval for MCMC samples.
#' @param ... If you have true values of the parameters, give it as "pars.true", then cp value return as 0 and 1 
#' for each parameter with respect to HPS interval includes true value or not.
#' 
#' @return A data frame containing the lower and upper bounds of the 95% HPD credible interval for each parameter. 
#' If true parameter values are provided, an additional column 'cp' indicates whether the true value falls within
#'  the interval (1) or not (0).
#'  
#' @examples
#' \dontrun{
#' # Assuming fit.MCMC.SMCM is fitted as "out", then
#' HPD.interval.mcmc(fit.mcmc.chains = out$fit, pars = "b.p", warmup = 1, thin = 1)
#' }
#' @export
HPD.interval.mcmc <- function(fit.mcmc.chains, pars, warmup, thin, ...  ){

  nchains                       <- length(fit.mcmc.chains)
  n_row                         <- nrow(fit.mcmc.chains[[1]][[pars]])
  n_col                         <- ncol(fit.mcmc.chains[[1]][[pars]])
  args                          <- list(...)
  pars.true                     <- if ("pars.true" %in% names(args)) args$pars.true else NA
  
  mcmc.chains                   <- array( unlist( lapply(fit.mcmc.chains, function(x) x[[pars]] ) ),  dim = c( n_row, n_col, nchains ) )
  
  dimnames(mcmc.chains)[[1]]    <- 1:n_row
  if(pars == "b.p"){
    dimnames(mcmc.chains)[[2]]  <- paste0(pars,"[", 0:(n_col-1),"]")
  }else{
    dimnames(mcmc.chains)[[2]]  <- paste0(pars,"[", 1:(n_col),"]")
  }
  dimnames(mcmc.chains)[[3]]    <- paste0("chain", 1:nchains) 
  
  mcmc.chains.list              <- mcmc.list(
    lapply(1:nchains, function(chain) { 
      df           <- as.data.frame(mcmc.chains[ , , chain])
      colnames(df) <- dimnames(mcmc.chains)[[2]]
      mcmc(df)
      # mcmc(as.data.frame( mcmc.chains[ , , chain ], col.names = dimnames(mcmc.chains)[[2]] ) )
    } ) )
  
  mcmc.chains.list.R1           <- window(mcmc.chains.list, start = warmup + 2, thin = thin)
  combined.mcmc.chains          <- do.call(rbind, mcmc.chains.list.R1) 
  
  hpd.interval                  <- HPDinterval( as.mcmc(combined.mcmc.chains), prob = 0.95 )
  hpd.interval                  <- as.data.frame(hpd.interval)           
  
  if( !is.na(pars.true[1]) ){
    for (i in 1:n_col) {
      hpd.interval[i, "cp"]     <- ifelse( (hpd.interval$lower[i]<= pars.true[i] & pars.true[i] <= hpd.interval$upper[i] ), 1, 0) } 
  }else{
    hpd.interval = hpd.interval
  }
  
  return(hpd.interval)
}