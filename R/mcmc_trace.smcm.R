#' @title Trace plot of MCMC draws for SMCM and SMCFM.
#' 
#' @description
#' Trace plot of posterior MCMC draws of the parameters including \eqn{\bm \beta}, \eqn{\bm b}, and 
#' \eqn{\bm \lambda} for SMCM and SMCFM based on our MCMC algorithms.
#' It is similar to \pkg{bayesplot} \link[bayesplot]{mcmc_trace} function but customized for our MCMC outputs. 
#'
#' @name mcmc_trace.smcm
#' 
#' @import ggplot2
#' 
#' @param fit.mcmc.chains A list of MCMC results for all the chains based on our MCMC algorithms.
#' @param pars Parameters can be "beta.p" or "b.p" or "lambda.p" or "theta.p"
#' @param warmup Number of warmup iterations.
#' @param thin Thinning interval for MCMC samples.
#' 
#'  
#' @return A trace plot of the specified parameters.
#' 
#' @examples
#' \dontrun{
#' # Assuming fit.MCMC.SMCM is fitted as "out", then
#' mcmc_trace.smcm(fit.mcmc.chains = out$fit, pars="b.p", warmup = 1, thin = 1 )
#' } 
#' @export
mcmc_trace.smcm <- function(fit.mcmc.chains, pars, warmup = 1, thin = 1 ){

  nchains     <- length(fit.mcmc.chains)
  n_row       <- nrow(fit.mcmc.chains[[1]][[pars]])
  n_col       <- ncol(fit.mcmc.chains[[1]][[pars]])
  
  mcmc.chains <-  array( unlist( lapply(fit.mcmc.chains, function(x) x[[pars]] ) ),  dim = c( n_row, n_col, nchains ) )
  
  dimnames(mcmc.chains)[[1]]    <- 1:n_row
  if(pars == "b.p"){
    dimnames(mcmc.chains)[[2]]  <-  paste0(pars,"[", 0:(n_col-1),"]")
  }else{
    dimnames(mcmc.chains)[[2]]  <-  paste0(pars,"[", 1:(n_col),"]")
  }
  dimnames(mcmc.chains)[[3]]    <- paste0("chain", 1:nchains) 
  
  mcmc.chains.list    <- coda::mcmc.list( lapply(1:nchains, function(chain) { 
    mcmc(as.data.frame( mcmc.chains[ , , chain], col.names = dimnames(mcmc.chains)[[2]] ) )
  } ) )
  mcmc.chains.list.R1 <- window(mcmc.chains.list, start = warmup + 2, thin = thin)
  
  bayesplot::color_scheme_set(c("#E66101" ,"#998EC3", "#542788", "#F1A340" ,"#D8DAEB", "#FEE0B6" ))
  # bayesplot::mcmc_trace(mcmc.chains.list, n_warmup = 0)
  trace_plot_smcm         <- bayesplot::mcmc_trace(mcmc.chains.list.R1, n_warmup = 0)
  
  trace_plot_smcm +
    scale_x_continuous(breaks = pretty(trace_plot_smcm$data$iteration), labels = pretty(seq(from=warmup+1, to=n_row, by=thin))[1:length(pretty(trace_plot_smcm$data$iteration))] ) 
  
}  
