// Rcpp versions of the expensive updates in SMCFM_fast_R1.R (mixture cure FRAILTY model)
//   UpdateRP_b_SMCFM_cpp      <->  UpdateRP.b.SMCFM.fast
//   UpdateRP_beta_SMCFM_cpp   <->  UpdateRP.beta.SMCFM.fast
//   UpdateBH_lambda_SMCFM_cpp <->  UpdateBH.lambda.SMCFM.fast
//   UpdateTheta_SMCFM_cpp     <->  UpdateTheta.SMCFM.fast
//
// The lasso hierarchy (UpdateTau / UpdateTauStar / UpdateSigma* / UpdateEta*) is shared with the SMCM
// model and is NOT repeated here (UpdateTau_cpp lives in SMCM_rcpp.cpp).
//
// Notation as in the R code:  H1 = H0 * exp(Xbeta),  l1p = log(1 + H1 / theta),
//   log S_frailty = -theta * l1p,  f_frailty = h0 exp(Xbeta) exp(-(theta + 1) * l1p).
//
// Random numbers are drawn in the same order as in the R code (rnorm, then runif only when
// logR is not NaN), so with the same seed the updates reproduce the R updates up to
// floating-point rounding.  `delta` is the 0/1 event indicator (R passes as.integer(event)).

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
using namespace Rcpp;

static inline double softplus_f(double x) {
  return std::max(x, 0.0) + std::log1p(std::exp(-std::fabs(x)));
}

static const double NEG_INF = -std::numeric_limits<double>::infinity();

// ---- log-likelihood pieces (same arithmetic and summation order as the R code) ----

// b-part:  -sum softplus(Zb) + sum_{event} Zb + sum_{cens} softplus(Zb + logS)
static double loglik_b_f(const double* zb, const double* logS, const int* d, int n) {
  long double sa = 0, sb = 0, sc = 0;
  for (int i = 0; i < n; ++i) {
    sa += softplus_f(zb[i]);
    if (d[i] == 1) sb += zb[i]; else sc += softplus_f(zb[i] + logS[i]);
  }
  double out = -static_cast<double>(sa);
  out += static_cast<double>(sb);
  out += static_cast<double>(sc);
  return out;
}

// beta-part:  sum_{event} [Xbeta - (theta+1) l1p] + sum_{cens} softplus(Zb - theta l1p)
static double loglik_beta_f(const double* xb, const double* zb, const double* H0,
                            double theta, const int* d, int n) {
  long double se = 0, sc = 0;
  for (int i = 0; i < n; ++i) {
    const double l1p = std::log1p(H0[i] * std::exp(xb[i]) / theta);
    if (d[i] == 1) se += (xb[i] - (theta + 1.0) * l1p); else sc += softplus_f(zb[i] - theta * l1p);
  }
  double out = 0.0;
  out += static_cast<double>(se);
  out += static_cast<double>(sc);
  return std::isfinite(out) ? out : NEG_INF;
}

// lambda-part:  sum_{event} [log h0 - (theta+1) l1p] + sum_{cens} softplus(Zb - theta l1p)
// `ex` = exp(Xbeta) does not change during the lambda update, so it is computed once by the caller.
static double loglik_lambda_f(const double* H0, const double* h0, const double* ex,
                              const double* zb, double theta, const int* d, int n) {
  for (int i = 0; i < n; ++i)
    if (!(h0[i] > 0.0) || !std::isfinite(h0[i])) return NEG_INF;
  long double se = 0, sc = 0;
  for (int i = 0; i < n; ++i) {
    const double l1p = std::log1p(H0[i] * ex[i] / theta);
    if (d[i] == 1) se += (std::log(h0[i]) - (theta + 1.0) * l1p);
    else           sc += softplus_f(zb[i] - theta * l1p);
  }
  double out = 0.0;
  out += static_cast<double>(se);
  out += static_cast<double>(sc);
  return std::isfinite(out) ? out : NEG_INF;
}

// theta-part (H1 fixed):  -(theta+1) sum_{event} l1p + sum_{cens} softplus(Zb - theta l1p)
static double loglik_theta_f(double theta, const double* H1, const double* zb, const int* d, int n) {
  long double se = 0, sc = 0;
  for (int i = 0; i < n; ++i) {
    const double l1p = std::log1p(H1[i] / theta);
    if (d[i] == 1) se += l1p; else sc += softplus_f(zb[i] - theta * l1p);
  }
  double out = -(theta + 1.0) * static_cast<double>(se);
  out += static_cast<double>(sc);
  return std::isfinite(out) ? out : NEG_INF;
}

// ---------------------------------------------------------------------------
// b update (incidence coefficients).  logS = -theta * log1p(H1 / theta) is computed in R.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List UpdateRP_b_SMCFM_cpp(const NumericMatrix& Z, const NumericVector& b, const NumericVector& Zb,
                          const NumericVector& logS, const IntegerVector& delta,
                          const NumericVector& sd_b, const NumericVector& proposal_sd) {
  const int n = Z.nrow();
  const int p = b.size();
  NumericVector bn = clone(b);
  NumericVector zb = clone(Zb);
  IntegerVector accept(p);
  std::vector<double> zbp(n);
  const double* ls = logS.begin();
  const int* d = delta.begin();

  double cur = loglik_b_f(zb.begin(), ls, d, n);

  for (int k = 0; k < p; ++k) {
    const double bprop = R::rnorm(bn[k], proposal_sd[k]);
    const double diff = bprop - bn[k];
    const double* zk = Z.begin() + static_cast<std::size_t>(k) * n;   // column k

    long double sa = 0, sb = 0, sc = 0;
    for (int i = 0; i < n; ++i) {
      const double v = zb[i] + zk[i] * diff;
      zbp[i] = v;
      sa += softplus_f(v);
      if (d[i] == 1) sb += v; else sc += softplus_f(v + ls[i]);
    }
    double prop = -static_cast<double>(sa);
    prop += static_cast<double>(sb);
    prop += static_cast<double>(sc);

    const double logR = prop - cur + R::dnorm(bprop, 0.0, sd_b[k], 1) -
                        R::dnorm(bn[k], 0.0, sd_b[k], 1);

    if (!std::isnan(logR) && std::log(R::runif(0.0, 1.0)) < logR) {
      bn[k] = bprop;
      std::copy(zbp.begin(), zbp.end(), zb.begin());
      cur = prop;
      accept[k] = 1;
    }
  }
  return List::create(_["b"] = bn, _["Zb"] = zb, _["loglik"] = cur, _["accept"] = accept);
}

// ---------------------------------------------------------------------------
// beta update (latency coefficients)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List UpdateRP_beta_SMCFM_cpp(const NumericMatrix& X, const NumericVector& beta, const NumericVector& Xbeta,
                             const NumericVector& Zb, const NumericVector& H0, double theta,
                             const IntegerVector& delta, const NumericVector& sd_beta,
                             const NumericVector& proposal_sd) {
  const int n = X.nrow();
  const int p = beta.size();
  NumericVector bn = clone(beta);
  NumericVector xb = clone(Xbeta);
  IntegerVector accept(p);
  std::vector<double> xbp(n);
  const double* zb = Zb.begin();
  const double* h0 = H0.begin();
  const int* d = delta.begin();

  double cur = loglik_beta_f(xb.begin(), zb, h0, theta, d, n);

  for (int k = 0; k < p; ++k) {
    const double bprop = R::rnorm(bn[k], proposal_sd[k]);
    const double diff = bprop - bn[k];
    const double* xk = X.begin() + static_cast<std::size_t>(k) * n;

    long double se = 0, sc = 0;
    for (int i = 0; i < n; ++i) {
      const double v = xb[i] + xk[i] * diff;
      xbp[i] = v;
      const double l1p = std::log1p(h0[i] * std::exp(v) / theta);
      if (d[i] == 1) se += (v - (theta + 1.0) * l1p); else sc += softplus_f(zb[i] - theta * l1p);
    }
    double prop = 0.0;
    prop += static_cast<double>(se);
    prop += static_cast<double>(sc);
    if (!std::isfinite(prop)) prop = NEG_INF;

    const double logR = prop - cur + R::dnorm(bprop, 0.0, sd_beta[k], 1) -
                        R::dnorm(bn[k], 0.0, sd_beta[k], 1);

    if (!std::isnan(logR) && std::log(R::runif(0.0, 1.0)) < logR) {
      bn[k] = bprop;
      std::copy(xbp.begin(), xbp.end(), xb.begin());
      cur = prop;
      accept[k] = 1;
    }
  }
  return List::create(_["beta"] = bn, _["Xbeta"] = xb, _["loglik"] = cur, _["accept"] = accept);
}

// ---------------------------------------------------------------------------
// lambda update (piecewise-constant baseline hazards), Gaussian random walk on log(lambda_k);
// the Jacobian log(lambda*) - log(lambda) is part of the ratio.  Prior lambda_k ~ Gamma(a, rate = b).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List UpdateBH_lambda_SMCFM_cpp(const NumericVector& lambda, const NumericVector& H0, const NumericVector& h0,
                               const NumericMatrix& A, const IntegerVector& interval_id,
                               const NumericVector& Xbeta, const NumericVector& Zb, double theta,
                               const IntegerVector& delta, double a, double b,
                               const NumericVector& proposal_sd) {
  const int n = A.nrow();
  const int J = lambda.size();
  NumericVector lam = clone(lambda);
  NumericVector H0c = clone(H0);
  NumericVector h0c = clone(h0);
  IntegerVector accept(J);
  std::vector<double> H0p(n), h0p(n), ex(n);
  const double* zb = Zb.begin();
  const int* d = delta.begin();
  const double scale = 1.0 / b;                       // R::dgamma uses scale = 1 / rate
  for (int i = 0; i < n; ++i) ex[i] = std::exp(Xbeta[i]);

  double cur = loglik_lambda_f(H0c.begin(), h0c.begin(), ex.data(), zb, theta, d, n);

  for (int k = 0; k < J; ++k) {
    const double lamc = lam[k];
    const double loglamc = std::log(lamc);
    const double loglamp = loglamc + R::rnorm(0.0, proposal_sd[k]);
    if (!std::isfinite(loglamp)) continue;
    const double lamp = std::exp(loglamp);
    if (!std::isfinite(lamp) || lamp <= 0.0) continue;

    const double diff = lamp - lamc;
    const double* Ak = A.begin() + static_cast<std::size_t>(k) * n;
    bool bad = false;
    for (int i = 0; i < n; ++i) {
      H0p[i] = H0c[i] + Ak[i] * diff;
      if (!std::isfinite(H0p[i]) || H0p[i] < 0.0) { bad = true; break; }
    }
    if (bad) continue;
    for (int i = 0; i < n; ++i) h0p[i] = (interval_id[i] == k + 1) ? lamp : h0c[i];

    const double prop = loglik_lambda_f(H0p.data(), h0p.data(), ex.data(), zb, theta, d, n);
    if (!std::isfinite(prop)) continue;

    const double lpc = R::dgamma(lamc, a, scale, 1);
    const double lpp = R::dgamma(lamp, a, scale, 1);
    const double logR = prop - cur + lpp - lpc + loglamp - loglamc;

    if (!std::isnan(logR) && std::log(R::runif(0.0, 1.0)) < logR) {
      lam[k] = lamp;
      std::copy(H0p.begin(), H0p.end(), H0c.begin());
      std::copy(h0p.begin(), h0p.end(), h0c.begin());
      cur = prop;
      accept[k] = 1;
    }
  }
  return List::create(_["lambda"] = lam, _["H0"] = H0c, _["h0"] = h0c,
                      _["loglik"] = cur, _["accept"] = accept);
}

// ---------------------------------------------------------------------------
// theta update (frailty parameter): Gaussian random walk on log(theta),
// prior theta ~ Gamma(shape = c_shape, rate = d_rate), H1 = H0 * exp(Xbeta) fixed.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List UpdateTheta_SMCFM_cpp(double theta, const NumericVector& H1, const NumericVector& Zb,
                           const IntegerVector& delta, double c_shape, double d_rate,
                           double proposal_sd) {
  const int n = H1.size();
  const int* d = delta.begin();
  const double cur = loglik_theta_f(theta, H1.begin(), Zb.begin(), d, n);

  const double logthp = std::log(theta) + R::rnorm(0.0, proposal_sd);
  const double thp = std::exp(logthp);
  if (!std::isfinite(thp) || thp <= 0.0) {
    return List::create(_["theta"] = theta, _["loglik"] = cur, _["accept"] = 0);
  }
  const double prop = loglik_theta_f(thp, H1.begin(), Zb.begin(), d, n);

  const double scale = 1.0 / d_rate;
  const double logR = prop - cur + R::dgamma(thp, c_shape, scale, 1) -
                      R::dgamma(theta, c_shape, scale, 1) + logthp - std::log(theta);

  if (!std::isnan(logR) && std::log(R::runif(0.0, 1.0)) < logR) {
    return List::create(_["theta"] = thp, _["loglik"] = prop, _["accept"] = 1);
  }
  return List::create(_["theta"] = theta, _["loglik"] = cur, _["accept"] = 0);
}
