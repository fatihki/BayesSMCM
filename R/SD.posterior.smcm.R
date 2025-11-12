#' @title Calculate posterior standard deviation of the estimates from SMCM or SMCFM fitting based on MCMC or Stan
#' 
#' @description 
#' This function calculates the posterior standard deviation of the parameter estimates from the fitted SMCM
#'  or SMCFM model based on MCMC or Stan.
#' It processes the MCMC chains or Stan outputs to compute the standard deviation for each parameter.
#' 
#' 
#' @name SD.posterior.smcm
#' 
#' @param out A list variable, an output from fit.MCMC.SMCM or fit.MCMC.SMCFM or any proposed Stan model.
#' @param stan.model Logical variable, if TRUE, it is from Stan, if FALSE, it is from our MCMC.
#' @param frailty Logical variable, if TRUE, it is cure frailty model, if FALSE, it is cure model.
#' 
#' @return A numeric vector containing the posterior standard deviations for each parameter.
#' 
#' @examples
#' \dontrun{
#' # Assuming fit.MCMC.SMCM is fitted as "out",  then
#' SD.posterior.smcm(out = out, stan.model = FALSE, frailty = FALSE)
#' }
#' 
#' @export
SD.posterior.smcm         <- function(out, stan.model = FALSE, frailty = FALSE){
  
  if(stan.model == FALSE){
    out$b.chains          <- do.call(rbind, out$b.chains)
    if( is.null(ncol(out$beta.chains[[1]])) == FALSE){
      out$beta.chains     <- do.call(rbind, out$beta.chains)
    }else{
      out$beta.chains     <- matrix(unlist(out$beta.chains), ncol = 1)
    }
    if( is.null(ncol(out$lambda.chains[[1]])) == FALSE ){
      out$lambda.chains   <- do.call(rbind, out$lambda.chains)
    }else{
      out$lambda.chains   <- matrix(unlist(out$lambda.chains), ncol = 1)
    }
    if(frailty == TRUE){
      out$theta.chains    <- unlist(out$theta.chains)
    }
  }
  
  b.chains.sd             <- apply(out$b.chains, 2, sd)
  beta.chains.sd          <- apply(out$beta.chains, 2, sd)
  lambda.chains.sd        <- apply(out$lambda.chains, 2, sd)
  out.sd                  <- c( SD.b.hat = b.chains.sd, SD.beta.hat = beta.chains.sd,  SD.lambda.hat = lambda.chains.sd )
  
  if(frailty == TRUE){
    theta.chains.sd  <- sd(out$theta.chains)
    out.sd                <- c( SD.b.hat = b.chains.sd, SD.beta.hat = beta.chains.sd, 
                                D.lambda.hat = lambda.chains.sd, SD.theta.hat = theta.chains.sd )
  }
  
  return(out.sd)
}
