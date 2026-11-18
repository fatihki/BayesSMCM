// Hierarchical semiparametric mixture cure FRAILTY model (HSMCFM)
//   Population survival: S_pop(t) = (1 - pi) + pi * S_u(t),   pi = logit^-1(Z b)
//   Uncured part: gamma frailty W ~ Gamma(theta, theta) (mean 1, variance 1/theta) acting on
//                 h(t | x, W) = W * lambda_j * exp(x' beta),  t in (s_j, s_{j+1}]
//   Marginal: S_u(t | x) = (1 + H0(t) exp(x' beta) / theta)^(-theta)
//             f_u(t | x) = h0(t) exp(x' beta) (1 + H0(t) exp(x' beta) / theta)^(-(theta + 1))
//   Priors: Bayesian-lasso (conditional Laplace) on b | sigmaSq and beta | sigmaStarSq,
//           sigmaSq, sigmaStarSq ~ Inverse-Gamma

functions {
  // Prior sd of b given the hierarchy. If the intercept (b[1]) is NOT penalised (as in the Rcpp samplers,
  // penalize_intercept = FALSE) it gets a N(0, intercept_sd^2) prior and only b[2:p2] use tauSq.
  vector b_sd(int p2, int penalize_intercept, real intercept_sd, vector tauSq, real sigmaSq) {
    vector[p2] sd_b;
    if (penalize_intercept == 1) {
      sd_b = sqrt(tauSq * sigmaSq);
    } else {
      sd_b[1] = intercept_sd;
      if (p2 > 1) sd_b[2:p2] = sqrt(tauSq * sigmaSq);
    }
    return sd_b;
  }

  // Per-subject log-likelihood of the mixture cure frailty model.
  //   delta  : event indicator (1 = event, 0 = censored)
  //   Zb     : cure-part linear predictor (logit of P(uncured))
  //   betaX  : survival-part linear predictor
  //   H0     : baseline cumulative hazard at t_i
  //   log_h0 : log baseline hazard at t_i
  //   theta  : gamma frailty parameter
  vector cure_frailty_loglik(vector delta, vector Zb, vector betaX, vector H0,
                             vector log_h0, real theta) {
    int n = rows(delta);
    vector[n] log_pz = log_inv_logit(Zb);                       // log(pi)
    vector[n] log_1mpz = log1m_inv_logit(Zb);                   // log(1 - pi)
    vector[n] log_term = log1p(H0 .* exp(betaX) / theta);       // log(1 + H0 e^{xb} / theta); exact at H0 = 0
    vector[n] ll;
    for (i in 1:n) {
      if (delta[i] == 1) {
        // event: pi * f_u(t)
        ll[i] = log_pz[i] + log_h0[i] + betaX[i] - (theta + 1) * log_term[i];
      } else {
        // censored: (1 - pi) + pi * S_u(t),  log S_u = -theta * log_term
        ll[i] = log_sum_exp(log_1mpz[i], log_pz[i] - theta * log_term[i]);
      }
    }
    return ll;
  }
}

data {
  int<lower=1> n;      // number of observations
  int<lower=1> p1;     // number of regressors for the survival part (X)
  int<lower=1> p2;     // number of regressors for the cure part (Z)
  int<lower=1> JJ;     // number of intervals
  vector[JJ+1] s;      // interval cut points: s[1] = 0 < s[2] < ... < s[JJ+1]
  vector[n] t;         // survival time
  vector[n] delta;     // status (1 = event, 0 = censored)
  matrix[n, p1] X;     // predictors, survival part
  matrix[n, p2] Z;     // predictors, cure part
  real<lower=0> lambda_a0;     // Gamma(shape, rate) prior for lambda_j
  real<lower=0> lambda_b0;
  real<lower=0> sigma_a0;      // Inverse-Gamma(shape, scale) prior for sigmaSq     (cure part)
  real<lower=0> sigma_b0;
  real<lower=0> sigmaStar_a0;  // Inverse-Gamma(shape, scale) prior for sigmaStarSq (survival part)
  real<lower=0> sigmaStar_b0;
  real<lower=0> r1;            // Gamma(shape, rate) prior for etaSq (in b)
  real<lower=0> delta1;
  real<lower=0> r2;            // Gamma(shape, rate) prior for etaStarSq (in beta)
  real<lower=0> delta2;
  real<lower=0> theta_a0;      // Gamma(shape, rate) prior for the frailty parameter theta
  real<lower=0> theta_b0;
  int<lower=0, upper=1> save_loglik;  // flag to save log-likelihood
  int<lower=0, upper=1> noncentered;  // 1 = non-centred parameterisation of b and beta, 0 = centred
  int<lower=0, upper=1> penalize_intercept;  // 0: b[1] (intercept) ~ N(0, intercept_sd^2), outside the lasso hierarchy
  real<lower=0> intercept_sd;
}

transformed data {
  // Same construction as in HSMCM_R1.stan / cumulative_hazard_matrix():
  //   expo[i, j] = (min(t_i, s_{j+1}) - s_j)_+  ->  H0(t_i) = (expo * lambda)[i]
  //   k[i] = interval containing t_i, intervals (s_j, s_{j+1}]; t_i = s[1] goes to interval 1.
  // All times must lie inside [s[1], s[JJ+1]].
  matrix[n, JJ] expo = rep_matrix(0, n, JJ);
  array[n] int k;
  int p2_pen = (penalize_intercept == 1) ? p2 : p2 - 1;   // number of penalised cure coefficients

  for (j in 1:JJ) {
    if (s[j+1] <= s[j])
      reject("interval bounds s must be strictly increasing; problem at j = ", j);
  }
  for (i in 1:n) {
    if (t[i] < s[1] || t[i] > s[JJ+1])
      reject("All observed times must lie inside the interval bounds; problem at i = ", i,
             ", t = ", t[i], ", range = [", s[1], ", ", s[JJ+1], "]");
    k[i] = 1;
    for (j in 1:JJ) {
      if (j > 1 && t[i] > s[j])
        k[i] = j;
      expo[i, j] = fmax(0, fmin(t[i], s[j+1]) - s[j]);
    }
  }
}

parameters {
  vector[p2] b_raw;                      // cure-part coefficients (centred) or standardised ones (non-centred)
  vector[p1] beta_raw;                   // survival-part coefficients (centred) or standardised ones (non-centred)
  vector<lower=0.00001>[JJ] lambda;      // piecewise constant baseline hazards
  real<lower=0> theta;                   // gamma frailty parameter
  real<lower=0> sigmaSq;
  real<lower=0> etaSq;
  vector<lower=0>[p2_pen] tauSq;             // auxiliary variances, conditional Laplace prior on b
  real<lower=0> sigmaStarSq;
  real<lower=0> etaStarSq;
  vector<lower=0>[p1] tauStarSq;         // auxiliary variances, conditional Laplace prior on beta
}

transformed parameters {
  // noncentered = 0: b, beta are the parameters themselves (centred; b_raw = b).
  // noncentered = 1: b = b_raw * sd, beta = beta_raw * sd, with b_raw, beta_raw ~ N(0, 1) a priori
  //                  (same prior on b and beta, better geometry when the data are weakly informative).
  vector[p2] b;
  vector[p1] beta;
  if (noncentered == 1) {
    b    = b_raw    .* b_sd(p2, penalize_intercept, intercept_sd, tauSq, sigmaSq);
    beta = beta_raw .* sqrt(tauStarSq * sigmaStarSq);
  } else {
    b    = b_raw;
    beta = beta_raw;
  }
}

model {
  // ---- priors ----
  lambda ~ gamma(lambda_a0, lambda_b0);
  theta  ~ gamma(theta_a0, theta_b0);

  // cure part
  sigmaSq ~ inv_gamma(sigma_a0, sigma_b0);
  etaSq   ~ gamma(r1, delta1);
  tauSq   ~ exponential(etaSq / 2);                 // rate = etaSq/2
  if (noncentered == 1) b_raw ~ std_normal();
  else                  b_raw ~ normal(0, b_sd(p2, penalize_intercept, intercept_sd, tauSq, sigmaSq));

  // survival part
  sigmaStarSq ~ inv_gamma(sigmaStar_a0, sigmaStar_b0);
  etaStarSq   ~ gamma(r2, delta2);
  tauStarSq   ~ exponential(etaStarSq / 2);
  if (noncentered == 1) beta_raw ~ std_normal();
  else                  beta_raw ~ normal(0, sqrt(tauStarSq * sigmaStarSq));

  // ---- likelihood ----
  {
    vector[n] H0 = expo * lambda;          // baseline cumulative hazard at t_i
    vector[n] log_h0 = log(lambda[k]);     // log baseline hazard at t_i
    target += sum(cure_frailty_loglik(delta, Z * b, X * beta, H0, log_h0, theta));
  }
}

generated quantities {
  vector[n] log_lik = rep_vector(0, n);
  if (save_loglik == 1) {
    log_lik = cure_frailty_loglik(delta, Z * b, X * beta, expo * lambda, log(lambda[k]), theta);
  }
}

