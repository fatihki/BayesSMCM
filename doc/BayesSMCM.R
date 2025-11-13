## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
#install.packages("remotes")
remotes::install_github("fatihki/BayesSMCM")

## -----------------------------------------------------------------------------
library(BayesSMCM)

## -----------------------------------------------------------------------------
# simulate data
n <- 300
b.true <- c( 0.4, 0.5, 0.1)
beta.true <- c(1, 0.2 )
baseline.hazard.rates <- 1
intervals <- c(0, 16)
prob.cov.X <- c(0.5)
prob.cov.Z <- c(0.5)
dat1 <- simSMCM(n, b.true, beta.true, baseline.hazard.rates, intervals, seed = 2026, cens.start = 30/365, cens.end = 32,  prob.cov.X, prob.cov.Z, same.cov = TRUE )
str(dat1)

## -----------------------------------------------------------------------------
# Bayesian settings:
nchains = 3; nIter   = 7500;   warmup = 2500; thin = 10
# priors for SMCMs
priorPar.smcm = list( r1 = 1, delta1 = 1e-04, r2 = 1, delta2 = 1e-04, a = 0.1, b = 0.1 )
# priors for SMCFMs
priorPar.smcfm = list( r1 = 1, delta1 = 1e-04, r2 = 1, delta2 = 1e-04, a = 0.1, b = 0.1, c = 0.1, d = 0.1 )

## -----------------------------------------------------------------------------
out.smcm.mcmc = fit.MCMC.SMCM ( data = dat1, hyperpar = priorPar.smcm, nchains, nIter, warmup, thin, mcmc.parallel = "parLapply",
                                standardize = FALSE, probs = 0, save_loglik = 1, seed = 2025 )

## -----------------------------------------------------------------------------
print_smcm( out.smcm.mcmc, stan.model = FALSE, frailty = FALSE )

## -----------------------------------------------------------------------------
out.smcfm.mcmc = fit.MCMC.SMCFM ( data = dat1, hyperpar = priorPar.smcfm, nchains, nIter, warmup, thin, mcmc.parallel = "parLapply",
                                  standardize = FALSE, probs = 0, save_loglik = 1, seed = 2025 )

## -----------------------------------------------------------------------------
print_smcm( out.smcfm.mcmc, stan.model = FALSE, frailty = TRUE )

## -----------------------------------------------------------------------------
priorPar.smcm.rstan = list( sigma_beta = 1000, sigma_b = 1000, a = 0.1, b = 0.1)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
out.smcm.rstan = fit.SMCM.RStan(data=dat1, hyperpar = priorPar.smcm.rstan, nchains, nIter, warmup, thin, standardize = FALSE,
                                probs = 0, save_loglik = 1, seed = 2025  )

## -----------------------------------------------------------------------------
print_smcm( out.smcm.rstan, stan.model = TRUE, frailty = FALSE)

## -----------------------------------------------------------------------------
priorPar.smcfm.rstan = list( sigma_beta = 1000, sigma_b = 1000,a = 0.1, b = 0.1,  c = 0.1, d = 0.1 )

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
out.smcfm.rstan =  fit.SMCFM.RStan(data=dat1, hyperpar =  priorPar.smcfm.rstan, nchains, nIter, warmup, thin, standardize = FALSE,
                                   probs = 0, save_loglik = 1, seed = 2025  )

## -----------------------------------------------------------------------------
print_smcm( out.smcfm.rstan, stan.model = TRUE, frailty = TRUE)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
out.hsmcm.rstan = fit.HSMCM.RStan(data=dat1, hyperpar =  priorPar.smcm , nchains, nIter, warmup, thin, standardize = FALSE,
                                  probs = 0, save_loglik = 1, seed = 2025  )

## -----------------------------------------------------------------------------
print_smcm( out.hsmcm.rstan, stan.model = TRUE, frailty = FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
out.hsmcfm.rstan = fit.HSMCFM.RStan(data=dat1, hyperpar = priorPar.smcfm, nchains, nIter, warmup, thin, standardize = FALSE,
                                    probs = 0, save_loglik = 1, seed = 2025  )

## -----------------------------------------------------------------------------
print_smcm( out.hsmcfm.rstan, stan.model = TRUE, frailty = TRUE)

## -----------------------------------------------------------------------------
data.E1690 = E1690
data.E1690.RFS <- data.E1690[data.E1690$failtime > 0, ]

real.data.RFS = list()
real.data.RFS$X =  as.matrix( cbind( trt=data.E1690.RFS$trt, 
                                     age=data.E1690.RFS$age, 
                                     sex=data.E1690.RFS$sex  ) )
head(real.data.RFS$X)
real.data.RFS$Z =  cbind(1, real.data.RFS$X)
head(real.data.RFS$Z)

real.data.RFS$observed_time = data.E1690.RFS$failtime 
real.data.RFS$delta         = data.E1690.RFS$rfscens 


## -----------------------------------------------------------------------------
# Bayesian settings:
nchains = 3; nIter   = 7500;   warmup = 2500; thin = 10
# priors for SMCMs
priorPar.smcm = list( r1 = 1, delta1 = 1e-04, r2 = 1, delta2 = 1e-04, a = 0.1, b = 0.1 )
# priors for SMCFMs
priorPar.smcfm = list( r1 = 1, delta1 = 1e-04, r2 = 1, delta2 = 1e-04, a = 0.1, b = 0.1, c = 0.1, d = 0.1 )

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
out.smcm.mcmc.E1690 = fit.MCMC.SMCM ( data = real.data.RFS, hyperpar = priorPar.smcm, nchains, nIter, warmup, thin, 
                                         mcmc.parallel = "parLapply", standardize = FALSE, probs = 0, save_loglik = 1, seed =  165251)

## -----------------------------------------------------------------------------
print_smcm( out.smcm.mcmc.E1690, stan.model = FALSE, frailty = FALSE )

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
out.hsmcm.rstan.E1690 = fit.HSMCM.RStan(data = real.data.RFS, hyperpar = priorPar.smcm , nchains, nIter, warmup, thin, 
                                        standardize = FALSE, probs = 0, save_loglik = 1, seed = 165251  )

## -----------------------------------------------------------------------------
print_smcm( out.hsmcm.rstan.E1690, stan.model = TRUE, frailty = FALSE )

## -----------------------------------------------------------------------------
# loading related packages
library(survival)
library(survminer)

## -----------------------------------------------------------------------------
# Preparing fitted model results for survival curves
b_mcmc_chains          <- lapply(out.smcm.mcmc.E1690$b.chains, coda::mcmc) 
beta_mcmc_chains       <- lapply(out.smcm.mcmc.E1690$beta.chains, coda::mcmc) 
lambda_mcmc_chains     <- lapply(out.smcm.mcmc.E1690$lambda.chains, coda::mcmc) 

b_mcmc_matrix          <- do.call(rbind, b_mcmc_chains)
beta_mcmc_matrix       <- do.call(rbind, beta_mcmc_chains)
lambda_mcmc_matrix     <- do.call(rbind, lambda_mcmc_chains)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
# --- KM estimate ---
fit.KM0    <- survfit(Surv(failtime, rfscens) ~ 1, data = data.E1690.RFS) 
km0_data   <- data.frame(
  time     = fit.KM0$time,
  surv     = fit.KM0$surv,
  n_event  = fit.KM0$n.event,
  n_censor = fit.KM0$n.censor
  )
km_censor0  <- subset(km0_data, n_censor > 0)

n           <- nrow(data.E1690.RFS)
S           <- nrow(b_mcmc_matrix)   # number of draws in the posterior sample
tgrid       <- fit.KM0$time
S_draws_avg <- matrix(NA, nrow = S, ncol = length(tgrid))

# function to scale covariate matrix which can include binary/categorical covariates, we only scale continuous covariates
scale.dummy.matrix <- function(X) {
  
  Xscaled = matrix(NA, nrow = nrow(X), ncol = ncol(X))
  
  for (j in 1:ncol(X)) {
    if ( all(X[,j] == floor(X[,j])) ) {
      Xscaled[,j] = X[,j]                                   # not scale if all the values are integer (for categorical variables)
    } else{
      Xscaled[,j] = scale(X[,j])
    }
  }
  colnames(Xscaled) <- colnames(X)
  return(Xscaled)
}

# scaled covariates
X           <- scale.dummy.matrix(real.data.RFS$X)
Z           <- scale.dummy.matrix(real.data.RFS$Z)

for (s in 1:S) {
  # Compute cumulative hazard and hazard of each tgrid[i]
  H0_h0           <- cumulative_hazard(tgrid, interval_bounds = out.smcm.mcmc.E1690$priors$s,
                                                  lambda = lambda_mcmc_matrix[s] ) 
  H0              <- H0_h0$H0 
  h0              <- H0_h0$h0 
  Zb              <- Z %*% b_mcmc_matrix[s,] 
  exp.Zb          <- exp(Zb)
  Xbeta           <- X %*% beta_mcmc_matrix[s,]
  exp.Xbeta       <- exp(Xbeta)
  H1              <- exp.Xbeta %*%  H0     # n x lengt(tgrid) matrix
  S_uncured       <- exp(-H1)              # n x T
  # Full cure model survival: n x T
  pi_z            <- exp.Zb/(1 + exp.Zb)
  S_i_t           <- sweep(S_uncured, 1, pi_z, FUN = "*")  # n x T
  S_i_t           <- sweep(S_i_t, 1, 1-pi_z, FUN = "+")    # n x T add pi_i to each row (broadcast row-wise)
  # Average over subjects (1 x T)
  S_draws_avg[s, ]<- colMeans(S_i_t)
}

# Posterior summaries across draws
S_mean_pop        <- apply(S_draws_avg, 2, mean)
S_low_pop         <- apply(S_draws_avg, 2, quantile, 0.025)
S_high_pop        <- apply(S_draws_avg, 2, quantile, 0.975)

model_data_pop    <- data.frame(
  time  = tgrid,
  surv  = S_mean_pop,
  lower = S_low_pop,
  upper = S_high_pop
)

## ----plot1, fig.cap="KM vs SMCM.MCMC Model."----------------------------------
real.data1.KM.model.plot <- ggplot() +
  # KM curve (step function)
  geom_step(data = km0_data, aes(x = time, y = surv),
            color = "black", linewidth = 0.8, linetype = "dashed") +
  # Add censor marks
  geom_point(data = km0_data,
             aes(x = time, y = surv),
             shape = 3,  # cross
             size = 2,
             stroke = 1) +
  # Model mean survival
  geom_line(data = model_data_pop,
            aes(x = time, y = surv),
            color = "#F8766D", linewidth = 1) +
  labs(x = "Time (years)", y = "Relapse-Free Survival Probability",
         title = "KM vs SMCM.MCMC Model") + 
  # Axis limits
  coord_cartesian(xlim = c(0, 7), ylim = c(0, 1)) +
  scale_x_continuous( breaks = seq(0, 7, by = 1)) +
  theme_classic(base_size = 14) 
real.data1.KM.model.plot

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
fit.KM <- survfit(Surv(failtime, rfscens) ~ trt, data = data.E1690.RFS)
km_data    <- data.frame(
  time     = fit.KM$time,
  surv     = fit.KM$surv,
  strata    = rep(names(fit.KM$strata), fit.KM$strata),  # repeat stratum label
  trt      = rep(c("0", "1"), fit.KM$strata),
  n_risk   = fit.KM$n.risk,
  n_event  = fit.KM$n.event,
  n_censor = fit.KM$n.censor
)
km_censor  <- subset(km_data, n_censor > 0)

tgrid <- sort(unique(fit.KM$time))

# For each draw, compute subject-specific survival (n x length(tgrid)), then average within trt group.
trt_levels              <- sort(unique(as.character(data.E1690.RFS$trt))) 
pop_draws_by_trt        <- list()
for (g in trt_levels){
  pop_draws_by_trt[[g]] <- matrix(NA, nrow = S, ncol = length(tgrid))
}

# Precompute group indices
idx_by_trt              <- lapply(trt_levels, function(x) which(as.character(data.E1690.RFS$trt) == x))
names(idx_by_trt)       <- trt_levels

S_draws_avg             <- matrix(NA, nrow = S, ncol = length(tgrid))
for (s in 1:S) {
  # Compute cumulative hazard and hazard of each tgrid[i]
  H0_h0            <- cumulative_hazard(tgrid, interval_bounds = out.smcm.mcmc.E1690$priors$s, lambda = lambda_mcmc_matrix[s] ) 
  H0               <- H0_h0$H0 
  h0               <- H0_h0$h0 
  Zb               <- Z %*% b_mcmc_matrix[s,] 
  exp.Zb           <- exp(Zb)
  Xbeta            <- X %*% beta_mcmc_matrix[s,]
  exp.Xbeta        <- exp(Xbeta)
  H1               <- exp.Xbeta %*%  H0    # n x lengt(tgrid) matrix
  S_uncured        <- exp(-H1)             # n x T
  # Full cure model survival: n x T
  pi_z             <- exp.Zb/(1 + exp.Zb)
  S_i_t            <- sweep(S_uncured, 1, pi_z, FUN = "*")  # n x T
  S_i_t            <- sweep(S_i_t, 1, 1-pi_z, FUN = "+")    # n x T add pi_i to each row (broadcast row-wise)
  # Average over subjects (1 x T)
  S_draws_avg[s, ] <- colMeans(S_i_t)
  
  # Alternative: average within each treatment group
  for (g in trt_levels) {
    idx <- idx_by_trt[[g]]
    if (length(idx) == 0) next
    pop_draws_by_trt[[g]][s, ] <- colMeans(S_i_t[idx, , drop = FALSE])
  }
}

# summarize across draws
model_pop_list        <- list()
for (g in trt_levels) {
  m                   <- pop_draws_by_trt[[g]]   # S x length(tgrid)
  mean_vec            <- apply(m, 2, mean)
  low_vec             <- apply(m, 2, quantile, 0.025)
  high_vec            <- apply(m, 2, quantile, 0.975)
  model_pop_list[[g]] <- data.frame(time = tgrid, surv = mean_vec, lower = low_vec, upper = high_vec, trt = g)
}
model_pop_df          <- dplyr::bind_rows(model_pop_list)


## ----message=FALSE, warning=FALSE---------------------------------------------
strata_cols   <- c("trt=0" = "black", "trt=1" = "gray")
trt_cols      <- c("0" = "#F8766D", "1" =  "#fcbbb6" ) #"#faa49e"

real.data1.KM.model.trt.plot <- ggplot() +
  # KM step lines by trt
  geom_step(data = km_data, aes(x = time, y = surv, color = strata, linetype = "KM"), linewidth = 0.8, show.legend = FALSE ) +
  # Add censor marks
  geom_point(data = km_censor,
             aes(x = time, y = surv, color = strata),
             shape = 3,  # cross
             size = 2,
             stroke = 1,  show.legend = FALSE ) +
  scale_color_manual(values = strata_cols,
                     name = "KM") +
  # --- New color scale for model curves ---
  ggnewscale::new_scale_color() +
  geom_line(data = model_pop_df, aes(x = time, y = surv, color = factor(trt),  linetype = "Model"), size = 1.2) +
  scale_color_manual(values = trt_cols,
                     labels = c("OBS", "High-dose IFN"),
                     name = "Treatment") +
  labs(x = "Time (years)", y = "Relapse-Free Survival Probability",
     title = "KM vs SMCM.MCMC Model by Treatment") +
  scale_linetype_manual(values = c("KM" = "dashed", "Model" = "solid")) +
  # Axis limits
  coord_cartesian(xlim = c(0, 7), ylim = c(0, 1)) +
  scale_x_continuous( breaks = seq(0, 7, by = 1)) +
  theme_classic(base_size = 14) +
  guides(linetype = "none") +
  theme(legend.position = "top")
real.data1.KM.model.trt.plot

