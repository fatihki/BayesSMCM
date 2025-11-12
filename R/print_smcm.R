#' @title Print the summary results from SMCM or SMCFM fitting based on MCMC or Stan
#' 
#' @description 
#' This function prints the summary results of the fitted SMCM or SMCFM model based on MCMC or Stan. 
#' It includes the posterior mean, posterior standard deviation, 95% HPD credible interval, and Rhat for each parameter.
#' The function handles both MCMC outputs from custom implementations and outputs from Stan models.
#' It also calculates the posterior standard deviation and HPD interval for the cure probability \eqn{\pi(z)}.
#' 
#' @name print_smcm
#' 
#' @importFrom coda mcmc as.mcmc mcmc.list gelman.diag HPDinterval
#' 
#' @param out A list variable, an output from fit.MCMC.SMCM or fit.MCMC.SMCFM or any proposed Stan model.
#' @param stan.model Logical variable, if TRUE, it is from Stan, if FALSE, it is from our MCMC.
#' @param frailty Logical variable, if TRUE, it is cure frailty model, if FALSE, it is cure model.
#' 
#' @return A data frame containing the summary results, including posterior mean, posterior SD, 
#' 95% HPD credible interval, and Rhat for all model parameters.
#'
#' @examples
#' \dontrun{
#' # Assuming fit.MCMC.SMCM is fitted as "out", then
#' print_smcm(out = out.smcm.mcmc, stan.model = FALSE, frailty = FALSE)
#' }
#'
#' @export
print_smcm <- function(out, stan.model, frailty){
  # if (is.null(out$mcmc.info)) {
  #   out$mcmc.info <- list(
  #     nchains = nchains,
  #     nIter = nIter,
  #     warmup = warmup,
  #     thin = thin
  #   )
  # }
  
  n_row         <- length(out$fit.results)-2
  df            <- array(NA, c(n_row,5) )
  length.lambda <- ifelse( (stan.model == FALSE), sum( grepl("lambda", names(out$fit.results)) ), ncol(out$lambda.chains) )
  
  if(frailty == FALSE) {
    rownames(df) <- c("b[0]", paste0("b[",1:ncol(out$data$X),"]"),  paste0("beta[",1:ncol(out$data$X),"]"), paste0("lambda[",1:length.lambda,"]") )
  }else{
    rownames(df) <- c("b[0]", paste0("b[",1:ncol(out$data$X),"]"),  paste0("beta[",1:ncol(out$data$X),"]"), paste0("lambda[",1:length.lambda,"]"),  "theta" )
  }
  
  colnames(df) <- c("Estimate(Mean)","SD","HPD lower","HPD upper", "Rhat")
  df[,1]       <- out$fit.results[1:n_row]
  df[,2]       <- SD.posterior.smcm(out, stan.model, frailty )
  
  if( out$mcmc.info$nchains > 1) {   
    df[,5]     <- out$Rhat[1:n_row] 
    }
  
  if (stan.model == FALSE) {
    cat( "\n-", ifelse( frailty == FALSE, "SMCM.MCMC", "SMCFM.MCMC"), "Model Results (J =", sum( grepl("lambda", names(out$fit.results)) ), ")" ,"\n" )
    hpd.lower.upper <- as.matrix( rbind(HPD.interval.mcmc(fit.mcmc.chains=out$fit, pars="b.p",      warmup=out$mcmc.info$warmup, thin=out$mcmc.info$thin ),
                                        HPD.interval.mcmc(fit.mcmc.chains=out$fit, pars="beta.p",   warmup=out$mcmc.info$warmup, thin=out$mcmc.info$thin ),
                                        HPD.interval.mcmc(fit.mcmc.chains=out$fit, pars="lambda.p", warmup=out$mcmc.info$warmup, thin=out$mcmc.info$thin) ), ncol=2)
    if(frailty == TRUE) {
      hpd.lower.upper <- as.matrix( rbind(hpd.lower.upper,
                                          HPD.interval.mcmc(fit.mcmc.chains=out$fit, pars="theta.p", warmup=out$mcmc.info$warmup, thin=out$mcmc.info$thin) ), ncol=2)
    }
    rownames(hpd.lower.upper) <- NULL
    colnames(hpd.lower.upper) <- NULL
    df[,3:4]    <-  hpd.lower.upper
    
  }else{
    cat( "\n-",  class(out) ,"Model Results (J =", sum( grepl("lambda", names(out$fit.results)) ), ")" ,"\n" )
    hpd.lower.upper <- as.matrix( rbind(HPD.credible.interval.stan.fit(out, pars = "b.p"),
                                        HPD.credible.interval.stan.fit(out, pars = "beta.p"),
                                        HPD.credible.interval.stan.fit(out, pars = "lambda.p")), ncol=2 )
    if(frailty == TRUE) {
      hpd.lower.upper <- as.matrix( rbind(hpd.lower.upper,
                                          HPD.credible.interval.stan.fit(out, pars="theta.p") ), ncol=2)
    }
    rownames(hpd.lower.upper) <- NULL
    colnames(hpd.lower.upper) <- NULL
    df[,3:4]    <-  hpd.lower.upper
  }
  
  
  if (stan.model == FALSE) {
    mcmc_chains_b     <- lapply( out$fit, function(x) { x$b.p[-c(1), ] } )      # removing the initial values!
    processed_chains_b<- lapply(mcmc_chains_b, function(chain) {
      # Apply burn-in and thinning
      thinned_chain           <- chain[(out$mcmc.info$warmup + 1):nrow(chain), ][seq(1, nrow(chain) - out$mcmc.info$warmup, by = out$mcmc.info$thin), ]
      # Create an mcmc object for each thinned chain
      mcmc(thinned_chain)
    })
    mcmc_list_b       <- mcmc.list(processed_chains_b)
    
    pz                <- lapply( mcmc_list_b, function (x) { 
      apply( x, 1, function (y){ mean(logit( scale.dummy.matrix(out$data$Z), y ))  } )
    } )
    mcmc_list_pz <- coda::mcmc.list(lapply(pz, mcmc))
    
    Rhat_pz           <- matrix( gelman.diag(mcmc_list_pz, autoburnin = F)$psrf[,1] , nrow=1 )
    sd_pz             <- mean(unlist(lapply(mcmc_list_pz, sd) ) )
    hpd.interval_pz   <- as.data.frame( HPDinterval( as.mcmc(unlist(mcmc_list_pz)), prob = 0.95 ) )
  }else{
    pz                <- apply(out$b.chains, 1, function (y){ mean(logit( scale.dummy.matrix(out$data$Z), y ))  } )
    sd_pz             <- sd(pz)
    hpd.interval_pz   <- as.data.frame(HPDinterval( coda::as.mcmc(pz), prob = 0.95 ) )
    pz_chains    <- mcmc.list(
      lapply(1:out$mcmc.info$nchains, function(i) {
        start         <- (i - 1) * (length(pz)/out$mcmc.info$nchains) + 1
        end           <- i * (length(pz)/out$mcmc.info$nchains)
        mcmc(as.matrix(pz)[start:end, ])
      })
    )
    Rhat_pz           <- matrix( gelman.diag(mcmc(pz_chains), autoburnin = F)$psrf[,1] , nrow=1 )
  }
  
  df                  <- rbind(df, c( out$fit.results["pz.hat"], sd_pz,  hpd.interval_pz$lower,  hpd.interval_pz$upper, Rhat_pz) )
  rownames(df)[nrow(df)] <- "pi(z)"
  print(round(df,4))
  return( invisible(df) )
}
