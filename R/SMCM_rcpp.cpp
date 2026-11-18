// Rcpp versions of the expensive updates in SMCM_fast_R1.R
//   UpdateRP_b_cpp      <->  UpdateRP.b.fast
//   UpdateRP_beta_cpp   <->  UpdateRP.beta.fast
//   UpdateBH_lambda_cpp <->  UpdateBH.lambda.fast
//   UpdateTau_cpp       <->  UpdateTau.fast / UpdateTauStar.fast
//
// The Metropolis updates draw random numbers in the same order as the R code
// (rnorm, then runif only when logR is not NaN), so with the same seed they
// should reproduce the R updates up to floating-point rounding.
// The inverse-Gaussian sampler is a different algorithm from SuppDists::rinvGauss,
// so tau draws agree in distribution, not draw by draw.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
using namespace Rcpp;

static inline double softplus_c(double x) {
  return std::max(x, 0.0) + std::log1p(std::exp(-std::fabs(x)));
}

// ---- log-likelihood pieces (same arithmetic and summation order as the R code) ----

static double loglik_b_c(const double* zb, const double* h1, const int* d, int n) {
  long double sa = 0, sb = 0, sc = 0;
  for (int i = 0; i < n; ++i) {
    sa += softplus_c(zb[i]);
    if (d[i] == 1) sb += zb[i]; else sc += softplus_c(zb[i] - h1[i]);
  }
  double out = -static_cast<double>(sa);
  out += static_cast<double>(sb);
  out += static_cast<double>(sc);
  return out;
}

static double loglik_beta_c(const double* xb, const double* zb, const double* H0,
                            const int* d, int n) {
  long double se = 0, sc = 0;
  for (int i = 0; i < n; ++i) {
    const double h1 = H0[i] * std::exp(xb[i]);
    if (d[i] == 1) se += (xb[i] - h1); else sc += softplus_c(zb[i] - h1);
  }
  double out = 0.0;
  out += static_cast<double>(se);
  out += static_cast<double>(sc);
  return std::isfinite(out) ? out : -std::numeric_limits<double>::infinity();
}

static double loglik_lambda_c(const double* H0, const double* h0, const double* xb,
                              const double* zb, const int* d, int n) {
  for (int i = 0; i < n; ++i)
    if (!(h0[i] > 0.0) || !std::isfinite(h0[i])) return -std::numeric_limits<double>::infinity();
  long double se = 0, sc = 0;
  for (int i = 0; i < n; ++i) {
    const double h1 = H0[i] * std::exp(xb[i]);
    if (d[i] == 1) se += (std::log(h0[i]) - h1); else sc += softplus_c(zb[i] - h1);
  }
  double out = 0.0;
  out += static_cast<double>(se);
  out += static_cast<double>(sc);
  return std::isfinite(out) ? out : -std::numeric_limits<double>::infinity();
}

// ---------------------------------------------------------------------------
// b update (incidence coefficients)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List UpdateRP_b_cpp(const NumericMatrix& Z, const NumericVector& b, const NumericVector& Zb,
                    const NumericVector& H1, const IntegerVector& delta,
                    const NumericVector& sd_b, const NumericVector& proposal_sd) {
  const int n = Z.nrow();
  const int p = b.size();
  NumericVector bn = clone(b);
  NumericVector zb = clone(Zb);
  IntegerVector accept(p);
  std::vector<double> zbp(n);
  const double* h1 = H1.begin();
  const int* d = delta.begin();

  double cur = loglik_b_c(zb.begin(), h1, d, n);

  for (int k = 0; k < p; ++k) {
    const double bprop = R::rnorm(bn[k], proposal_sd[k]);
    const double diff = bprop - bn[k];
    const double* zk = Z.begin() + static_cast<std::size_t>(k) * n;   // column k

    long double sa = 0, sb = 0, sc = 0;
    for (int i = 0; i < n; ++i) {
      const double v = zb[i] + zk[i] * diff;
      zbp[i] = v;
      sa += softplus_c(v);
      if (d[i] == 1) sb += v; else sc += softplus_c(v - h1[i]);
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
List UpdateRP_beta_cpp(const NumericMatrix& X, const NumericVector& beta, const NumericVector& Xbeta,
                       const NumericVector& Zb, const NumericVector& H0, const IntegerVector& delta,
                       const NumericVector& sd_beta, const NumericVector& proposal_sd) {
  const int n = X.nrow();
  const int p = beta.size();
  NumericVector bn = clone(beta);
  NumericVector xb = clone(Xbeta);
  IntegerVector accept(p);
  NumericVector H1(n);
  std::vector<double> xbp(n), h1p(n);
  const double* zb = Zb.begin();
  const double* h0 = H0.begin();
  const int* d = delta.begin();

  double cur = loglik_beta_c(xb.begin(), zb, h0, d, n);
  for (int i = 0; i < n; ++i) H1[i] = h0[i] * std::exp(xb[i]);

  for (int k = 0; k < p; ++k) {
    const double bprop = R::rnorm(bn[k], proposal_sd[k]);
    const double diff = bprop - bn[k];
    const double* xk = X.begin() + static_cast<std::size_t>(k) * n;

    long double se = 0, sc = 0;
    for (int i = 0; i < n; ++i) {
      const double v = xb[i] + xk[i] * diff;
      const double h = h0[i] * std::exp(v);
      xbp[i] = v;
      h1p[i] = h;
      if (d[i] == 1) se += (v - h); else sc += softplus_c(zb[i] - h);
    }
    double prop = 0.0;
    prop += static_cast<double>(se);
    prop += static_cast<double>(sc);
    if (!std::isfinite(prop)) prop = -std::numeric_limits<double>::infinity();

    const double logR = prop - cur + R::dnorm(bprop, 0.0, sd_beta[k], 1) -
                        R::dnorm(bn[k], 0.0, sd_beta[k], 1);

    if (!std::isnan(logR) && std::log(R::runif(0.0, 1.0)) < logR) {
      bn[k] = bprop;
      std::copy(xbp.begin(), xbp.end(), xb.begin());
      std::copy(h1p.begin(), h1p.end(), H1.begin());
      cur = prop;
      accept[k] = 1;
    }
  }
  return List::create(_["beta"] = bn, _["Xbeta"] = xb, _["H1"] = H1,
                      _["loglik"] = cur, _["accept"] = accept);
}

// ---------------------------------------------------------------------------
// lambda update (piecewise-constant baseline hazards), log-scale random walk
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List UpdateBH_lambda_cpp(const NumericVector& lambda, const NumericVector& H0, const NumericVector& h0,
                         const NumericMatrix& A, const IntegerVector& interval_id,
                         const NumericVector& Xbeta, const NumericVector& Zb, const IntegerVector& delta,
                         double a, double b, const NumericVector& proposal_sd) {
  const int n = A.nrow();
  const int J = lambda.size();
  NumericVector lam = clone(lambda);
  NumericVector H0c = clone(H0);
  NumericVector h0c = clone(h0);
  IntegerVector accept(J);
  std::vector<double> H0p(n), h0p(n);
  const double* xb = Xbeta.begin();
  const double* zb = Zb.begin();
  const int* d = delta.begin();
  const double scale = 1.0 / b;                       // R::dgamma uses scale = 1 / rate

  double cur = loglik_lambda_c(H0c.begin(), h0c.begin(), xb, zb, d, n);

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

    const double prop = loglik_lambda_c(H0p.data(), h0p.data(), xb, zb, d, n);
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
// Bayesian-lasso scale mixture: tau^2_k = 1 / gamma_k,
// gamma_k ~ InverseGaussian(mean = sqrt(eta2 * sigma2) / |b_k|, shape = eta2)
// Michael-Schucany-Haas sampler, written in the numerically stable form.
// ---------------------------------------------------------------------------
static double rinvgauss_c(double mu, double lambda) {
  const double z = R::rnorm(0.0, 1.0);
  const double y = z * z;
  const double t = mu * y / (2.0 * lambda);
  const double x = mu / (1.0 + t + std::sqrt(t * t + 2.0 * t));
  const double u = R::runif(0.0, 1.0);
  return (u <= mu / (mu + x)) ? x : mu * mu / x;
}

// [[Rcpp::export]]
NumericVector UpdateTau_cpp(const NumericVector& b, double sigmaSq, double etaSq, double eps = 1e-12) {
  const int p = b.size();
  NumericVector out(p);
  const double root = std::sqrt(etaSq * sigmaSq);
  for (int k = 0; k < p; ++k) {
    const double nu = root / std::max(std::fabs(b[k]), eps);
    double g = rinvgauss_c(nu, etaSq);
    for (int tries = 0; (!std::isfinite(g) || g <= 0.0) && tries < 100; ++tries)
      g = rinvgauss_c(nu, etaSq);
    out[k] = 1.0 / g;
  }
  return out;
}
