//Updated hierarchical semiparametric.mixture.cure.frailty.model
data {
  int<lower=1> n; //number of observations
  int<lower=1> p1; //number of regressors for the survival part (X)
  int<lower=1> p2; //number of regressors for the cure part (Z)
  int<lower=1> JJ; //number of intervals
  vector[JJ+1] s; //interval cut points
  vector[n] t; // survival time
  vector[n] delta; //status
  matrix[n, p1] X; //predictor related to the survival part
  matrix[n, p2] Z; //predictor related to the cure part
  real lambda_a0; // parameters of lambda prior
  real lambda_b0;
  real<lower=0> r1;  // Gamma parameter for eta^2 in b
  real<lower=0> delta1;  // Gamma parameter for eta^2 in b
  real<lower=0> r2;  // Gamma parameter for eta_star^2 in beta
  real<lower=0> delta2;  // Gamma parameter for eta_star^2 in beta
  real theta_a0; // parameters of Gamma frailty prior
  real theta_b0;
  int<lower=0, upper=1> save_loglik; // Flag to save log-likelihood
}
// The parameters accepted by the model. Our model
parameters {
  vector[p2] b; // regression coefficients for the cure rate
  vector[p1] beta; // regression coefficients for the survival part
  vector<lower=0.00001>[JJ] lambda; //  piecewise constant baseline hazards
  real<lower=0> theta; // Gamma frailty parameter
  real<lower=0> sigmaSq;           // Error term standard deviation
  real<lower=0> etaSq;           
  vector<lower=0>[p2] tauSq;        // Auxiliary variable for conditional Laplace prior
  real<lower=0> sigmaStarSq;           // Error term standard deviation
  real<lower=0> etaStarSq;  
  vector<lower=0>[p1] tauStarSq;        // Auxiliary variable for conditional Laplace prior
}
// The model to be estimated. We model the output
model {
  real H0;
  int current_interval;
  real log_term;
  real log_uncured_survival;
  vector[n] pz = inv_logit(Z * b);
  vector[n] betaX = X * beta;
  vector[n] exp_betaX = exp(betaX);
   // priors
    lambda ~ gamma( lambda_a0, lambda_b0);
    sigmaSq  ~ uniform(0, 1000);                     // Approximates non-informative prior on 1/sigma^2 or use uniform(0, 1000) cauchy(0, 10); 
    etaSq  ~ gamma(r1, delta1);
      // Conditional Laplace prior on b given sigma
  for (j in 1:p2) {
    tauSq[j] ~ exponential(etaSq/2);              // Exponential prior for auxiliary variable
    b[j]   ~ normal(0, sqrt(tauSq[j] * sigmaSq) );   // Equivalent to conditional Laplace prior on beta | sigma_star^2 
  }
    sigmaStarSq  ~ uniform(0, 1000);                     // Approximates non-informative prior on 1/sigma_star^2  or use uniform(0, 1000) cauchy(0, 10); 
    etaStarSq  ~ gamma(r2, delta2);
      // Conditional Laplace prior on beta given sigmaStar
  for (j in 1:p1) {
    tauStarSq[j] ~ exponential(etaStarSq/2);              // Exponential prior for auxiliary variable
    beta[j]    ~ normal(0, sqrt(tauStarSq[j] * sigmaStarSq) );   // Equivalent to conditional Laplace prior on beta | sigma_star^2 
  }
   theta ~  gamma( theta_a0, theta_b0); //uniform(0, 1000); //  prior for the Gamma frailty parameter  uniform(0, 1000); //
     // likelihood
  for (i in 1:n) {
     // Compute cumulative hazard
    H0 = 0;
    current_interval = 0;
    for (j in 1:JJ) {
      if (t[i] > s[j] && t[i] <= s[j+1]) {
         // Add the hazard contribution for the current interval
        H0 += lambda[j] * (t[i] - s[j]);
        current_interval = j;
        break;
      } else if (t[i] > s[j+1]) {
        // Add the full hazard contribution for earlier intervals
        H0 += lambda[j] * (s[j+1] - s[j]);
      }
    }
    // real exp_beta_X = exp(dot_product(beta, X[i]));
    log_term = log1p(H0 * exp_betaX[i] / theta + 1e-8); // Stabilized log-term
    log_uncured_survival = -theta * log_term;
    if (delta[i] == 1) { // Event observed
      // Use the hazard rate for the correct interval
      target += log(pz[i]) + log(lambda[current_interval]) + betaX[i] - (theta + 1) * log_term;
    } else { // Censored
      target += log_sum_exp( log(1-pz[i]), log(pz[i]) + log_uncured_survival ) ; // log_sum_exp(a,b) = log(exp(a) + exp(b) )
    }
  }
}
// the computed log_lik matrix
generated quantities {
    vector[n] log_lik;
      if (save_loglik == 1) {
    vector[n] pz = inv_logit(Z * b);
    vector[n] betaX = X * beta;
    vector[n] H0;
    vector[n] h0;
    vector[n] log_term;
    int current_interval;
        for (i in 1:n) {
    // Compute cumulative hazard
    H0[i] = 0;
    for (j in 1:JJ) {
      if (t[i] > s[j] && t[i] <= s[j+1]) {
         // Add the hazard contribution for the current interval
        H0[i] = H0[i] + lambda[j] * (t[i] - s[j]);
        current_interval  = j;
        break;
      } else if (t[i] > s[j+1]) {
        // Add the full hazard contribution for earlier intervals
        H0[i] = H0[i] + lambda[j] * (s[j+1] - s[j]);
      }
    }
    h0[i] = lambda[current_interval];
    log_term[i] = log1p(H0[i] * exp(betaX)[i] / theta + 1e-8); // Stabilized log-term
        if (delta[i] == 1) { // Event observed
          // Use the hazard rate for the correct interval
          log_lik[i] = log(pz[i]) + log(h0[i]) + betaX[i] - (theta + 1) * log_term[i];
        } else { // Censored
          log_lik[i] = log_sum_exp( log(1-pz[i]), log(pz[i]) + (-theta * log_term[i]) ) ; // log_sum_exp(a,b) = log(exp(a) + exp(b) )
          }
        }
    } else {
        log_lik = rep_vector(0, n); // If log-lik is not needed, return zeros
    }
}
