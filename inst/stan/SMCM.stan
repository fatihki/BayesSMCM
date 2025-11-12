// updated semiparametric.mixture.cure.model
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
  real beta_sigma; // variance of beta
  real b_sigma; // variance of b
  int<lower=0, upper=1> save_loglik; // Flag to save log-likelihood
}
// The parameters accepted by the model. Our model
parameters {
  vector[p2] b; // regression coefficients for the cure rate
  vector[p1] beta; // regression coefficients for the survival part
  vector<lower=0.00001>[JJ] lambda; //  piecewise constant baseline hazards
}
// The model to be estimated. We model the output
model {
  real H0;
  int current_interval;
  real uncured_survival;
  vector[n] pz = inv_logit(Z * b);
  vector[n] betaX = X * beta;
  vector[n] exp_betaX = exp(betaX);
   // priors
     lambda ~ gamma( lambda_a0, lambda_b0);
     beta ~ normal( 0, beta_sigma);
     b ~ normal( 0, b_sigma);
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
    uncured_survival = exp(-H0 * exp_betaX[i]);
    if (delta[i] == 1) { // Event observed
      // Use the hazard rate for the correct interval
      target += log(pz[i]) + log(lambda[current_interval]) + betaX[i] + log(uncured_survival);
    } else { // Censored
      target += log_sum_exp( log(1-pz[i]), log(pz[i]) + log(uncured_survival) ) ; // log_sum_exp(a,b) = log(exp(a) + exp(b) )
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
    vector[n] uncured_survival;
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
    uncured_survival[i] = exp(-H0[i] * exp(betaX)[i]);
        if (delta[i] == 1) { // Event observed
          // Use the hazard rate for the correct interval
          log_lik[i] = log(pz[i]) + log(h0[i]) + betaX[i] + log(uncured_survival[i]);
        } else { // Censored
          log_lik[i] = log_sum_exp( log(1-pz[i]), log(pz[i]) + log(uncured_survival[i]) ) ; // log_sum_exp(a,b) = log(exp(a) + exp(b) )
          }
        }
    } else {
        log_lik = rep_vector(0, n); // If log-lik is not needed, return zeros
    }
}
