#' @title Gelman and Rubin's convergence diagnostic for MCMC chains of SMCM and SMCFM.
#' 
#' @description
#' The ‘potential scale reduction factor’ is calculated by using \pkg{coda} via \link[coda]{gelman.diag}.
#' 
#' @name Rhat.mcmc
#' 
#' @importFrom coda mcmc as.mcmc mcmc.list gelman.diag
#' 
#' @param fit A list variable \code{fit} from the fitted model of MCMC.SMCM or MCMC.SMCFM output which includes posterior samples of
#'  all the parameters over chains.
#' @param warmup Number of warmup iterations.
#' @param thin Thinning interval for MCMC samples.
#' 
#'  
#' @return Rhat values for all model parameters \eqn{\bm b}, \eqn{\bm \beta}, \eqn{\bm \lambda} and/or \eqn{\theta}.
#' 
#' @examples
#' \dontrun{
#' # Assuming fit.MCMC.SMCM or fit.MCMC.SMCFM is fitted as "out", then  
#' Rhat_values <- Rhat.mcmc(fit = out$fit, warmup, thin)
#' print(Rhat_values)
#' } 
#' 
#' @export
Rhat.mcmc                   <- function(fit, warmup = 1, thin = 1 ){
  
  # for b  
  mcmc_chains_b             <- lapply( fit, function(x) { x$b.p[-c(1), ] } )    # removing the initial values!
  processed_chains_b        <- lapply(mcmc_chains_b, function(chain) {
    # Apply burn-in and thinning
    thinned_chain           <- chain[(warmup + 1):nrow(chain), ][seq(1, nrow(chain) - warmup, by = thin), ]
    # Create an mcmc object for each thinned chain
    mcmc(thinned_chain)
  })
  
  mcmc_list_b               <- mcmc.list(processed_chains_b)
  Rhat_b                    <- matrix( gelman.diag(mcmc_list_b, autoburnin = F)$psrf[,1] , nrow=1 )
  colnames(Rhat_b )         <- paste0("b[", 1:length(Rhat_b), "]" ) 
  
  # for beta
  mcmc_chains_beta           <- lapply( fit, function(x) { x$beta.p[-c(1), ] } ) # removing the initial values!
  if( ncol(fit[[1]]$beta.p) == 1 ){
    processed_chains_beta    <- lapply(mcmc_chains_beta, function(chain) {
      thinned_chain          <- chain[(warmup + 1):length(chain) ][seq(1, length(chain) - warmup, by = thin) ]
      mcmc(thinned_chain)
    })
  }else{
    processed_chains_beta    <- lapply(mcmc_chains_beta, function(chain) {
      thinned_chain          <- chain[(warmup + 1):nrow(chain), ][seq(1, nrow(chain) - warmup, by = thin), ]
      mcmc(thinned_chain)
    })
  }
  
  mcmc_list_beta             <- mcmc.list(processed_chains_beta)
  Rhat_beta                  <- matrix( gelman.diag(mcmc_list_beta, autoburnin = F)$psrf[,1] , nrow=1 )
  colnames(Rhat_beta )       <- paste0("beta[", 1:length(Rhat_beta), "]" ) 
  
  # for lambda
  if( ncol(fit[[1]]$lambda.p) == 1 ){
    mcmc_chains_lambda       <- lapply( fit, function(x) { x$lambda.p[-c(1), ] } ) # removing the initial values!
    processed_chains_lambda  <- lapply(mcmc_chains_lambda, function(chain) {
      thinned_chain          <- chain[(warmup + 1):length(chain) ][seq(1, length(chain) - warmup, by = thin) ]
      mcmc(thinned_chain)
    })
  }else{
    mcmc_chains_lambda       <- lapply( fit, function(x) { x$lambda.p[-c(1), ] } ) # removing the initial values!
    processed_chains_lambda  <- lapply(mcmc_chains_lambda, function(chain) {
      thinned_chain          <- chain[(warmup + 1):nrow(chain), ][seq(1, nrow(chain) - warmup, by = thin), ]
      mcmc(thinned_chain)
    })
  }
  
  mcmc_list_lambda           <- mcmc.list(processed_chains_lambda)
  Rhat_lambda                <- matrix( gelman.diag(mcmc_list_lambda, autoburnin = F)$psrf[,1] , nrow=1 )
  colnames(Rhat_lambda )     <- paste0("lambda[", 1:length(Rhat_lambda), "]" ) 
  
  Rhat.all                   <- cbind(Rhat_b, Rhat_beta, Rhat_lambda) 
  
  # for theta (frailty parameter)
  if(!is.null(fit[[1]]$theta.p)){
    mcmc_chains_theta         <- lapply( fit, function(x) { x$theta.p[-c(1), ] } )    # removing the initial values!
    processed_chains_theta    <- lapply(mcmc_chains_theta, function(chain) {
      thinned_chain           <- chain[(warmup + 1):length(chain) ][seq(1, length(chain) - warmup, by = thin) ]
      mcmc(thinned_chain)
    })
    
    mcmc_list_theta           <- mcmc.list(processed_chains_theta)
    Rhat_theta                <- matrix( gelman.diag(mcmc_list_theta, autoburnin = F)$psrf[,1] , nrow=1 )
    colnames(Rhat_theta )     <- "theta"
    
    Rhat.all                  <- cbind(Rhat_b, Rhat_beta, Rhat_lambda, Rhat_theta) 
  }
  
  return( Rhat.all )
}