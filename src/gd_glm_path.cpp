/* Compiled Stage A path of the univariate GLM fitter in GD path mode.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 */

#include <Rcpp.h>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace Rcpp;

NumericVector Knotnew(NumericVector weights, NumericVector residuals, NumericVector x,
                      NumericVector dcum, NumericVector oldknots, double tol,
                      int support_order);

namespace {

// Family and link codes; see .gd_glm_fast_codes() in R/UnivariateFitter.R.
enum Family { GAUSSIAN = 0, POISSON = 1, BINOMIAL = 2, GAMMA = 3, INVERSE_GAUSSIAN = 4 };
enum Link { IDENTITY = 0, LOG = 1, LOGIT = 2, INVERSE = 3, INVERSE_SQUARE = 4 };

const double THRESH = 30.0, MTHRESH = -30.0, INVEPS = 1.0 / DBL_EPSILON;

// The link functions of stats::make.link(), including its clamping.
double linkfun(int link, double mu) {
  switch (link) {
  case IDENTITY: return mu;
  case LOG: return std::log(mu);
  case LOGIT: return std::log(mu / (1.0 - mu));
  case INVERSE: return 1.0 / mu;
  default: return 1.0 / (mu * mu);
  }
}
double linkinv(int link, double eta) {
  switch (link) {
  case IDENTITY: return eta;
  case LOG: return std::max(std::exp(eta), DBL_EPSILON);
  case LOGIT: {
    const double tmp = eta < MTHRESH ? DBL_EPSILON : (eta > THRESH ? INVEPS : std::exp(eta));
    return tmp / (1.0 + tmp);
  }
  case INVERSE: return 1.0 / eta;
  default: return 1.0 / std::sqrt(eta);
  }
}
double mu_eta(int link, double eta) {
  switch (link) {
  case IDENTITY: return 1.0;
  case LOG: return std::max(std::exp(eta), DBL_EPSILON);
  case LOGIT: {
    const double opexp = 1.0 + std::exp(eta);
    return (eta > THRESH || eta < MTHRESH) ? DBL_EPSILON : std::exp(eta) / (opexp * opexp);
  }
  case INVERSE: return -1.0 / (eta * eta);
  default: return -1.0 / (2.0 * std::pow(eta, 1.5));
  }
}
bool valideta(int link, const std::vector<double>& eta) {
  if (link == INVERSE) {
    for (double e : eta) if (!R_FINITE(e) || e == 0.0) return false;
  } else if (link == INVERSE_SQUARE) {
    for (double e : eta) if (!R_FINITE(e) || !(e > 0.0)) return false;
  }
  return true;
}
double variance(int family, double mu) {
  switch (family) {
  case GAUSSIAN: return 1.0;
  case POISSON: return mu;
  case BINOMIAL: return mu * (1.0 - mu);
  case GAMMA: return mu * mu;
  default: return mu * mu * mu;
  }
}
bool validmu(int family, const std::vector<double>& mu) {
  if (family == POISSON || family == GAMMA) {
    for (double m : mu) if (!R_FINITE(m) || !(m > 0.0)) return false;
  } else if (family == BINOMIAL) {
    for (double m : mu) if (!R_FINITE(m) || !(m > 0.0 && m < 1.0)) return false;
  }
  return true;
}
double y_log_y(double y, double mu) { return y != 0.0 ? y * std::log(y / mu) : 0.0; }
// Unit deviances of the stats families, summed as R's long-double sum().
double deviance(int family, const double* y, const std::vector<double>& mu,
                const double* wt, int n) {
  long double total = 0.0;
  for (int i = 0; i < n; ++i) {
    const double yi = y[i], m = mu[i], w = wt[i];
    double r;
    switch (family) {
    case GAUSSIAN: r = w * ((yi - m) * (yi - m)); break;
    case POISSON: r = 2.0 * (yi > 0.0 ? w * (yi * std::log(yi / m) - (yi - m)) : m * w); break;
    case BINOMIAL: r = 2.0 * w * (y_log_y(yi, m) + y_log_y(1.0 - yi, 1.0 - m)); break;
    case GAMMA: r = -2.0 * w * (std::log(yi == 0.0 ? 1.0 : yi / m) - (yi - m) / m); break;
    default: r = w * ((yi - m) * (yi - m)) / (yi * (m * m)); break;
    }
    total += r;
  }
  return (double) total;
}

// Linear-spline (hat-function) basis on the breakpoints t: row i is
// (1 - u[i]) at column cell[i] and u[i] at cell[i] + 1, or zero outside [t0, tm],
// as splineDesign(ord = 2, outer.ok = TRUE).
void hat_basis(const double* x, int n, const std::vector<double>& t,
               std::vector<int>& cell, std::vector<double>& u) {
  const int m = t.size();
  cell.assign(n, -1); u.assign(n, 0.0);
  for (int i = 0; i < n; ++i) {
    const double xi = x[i];
    if (xi < t[0] || xi > t[m - 1]) continue;
    int k = (int) (std::upper_bound(t.begin(), t.end(), xi) - t.begin()) - 1;
    if (k > m - 2) k = m - 2;
    cell[i] = k;
    u[i] = (xi - t[k]) / (t[k + 1] - t[k]);
  }
}

// Weighted least squares for the hat basis: tridiagonal normal equations,
// solved by an LDL' sweep. Declines singular or poorly conditioned systems.
bool hat_wls(const std::vector<int>& cell, const std::vector<double>& u,
             const std::vector<double>& z, const std::vector<double>& w, int m,
             std::vector<double>& theta) {
  std::vector<double> d(m, 0.0), e(m - 1, 0.0), b(m, 0.0);
  const int n = cell.size();
  for (int i = 0; i < n; ++i) {
    const int k = cell[i];
    const double wi = w[i];
    if (k < 0 || wi == 0.0) continue;
    const double ui = u[i], ai = 1.0 - ui;
    d[k] += wi * ai * ai;
    d[k + 1] += wi * ui * ui;
    e[k] += wi * ai * ui;
    b[k] += wi * ai * z[i];
    b[k + 1] += wi * ui * z[i];
  }
  const double rel_tol = 1e-10;
  if (!(d[0] > 0.0)) return false;
  for (int k = 1; k < m; ++k) {
    const double l = e[k - 1] / d[k - 1];
    const double next = d[k] - l * e[k - 1];
    if (!(d[k] > 0.0) || !(next > rel_tol * d[k])) return false;
    b[k] -= l * b[k - 1];
    d[k] = next;
  }
  theta.assign(m, 0.0);
  theta[m - 1] = b[m - 1] / d[m - 1];
  for (int k = m - 2; k >= 0; --k) theta[k] = (b[k] - e[k] * theta[k + 1]) / d[k];
  for (double v : theta) if (!R_FINITE(v)) return false;
  return true;
}

struct Fit {
  std::vector<double> theta, residuals, weights;
  double deviance;
};

// Outcome of one glm.fit() run: a fit (converged or not), an error that R
// would raise, or a branch not reproduced here (handed back to R).
enum Outcome { FIT, ERROR, UNSUPPORTED };

// glm.fit(x, y, weights, mustart = mustart, offset) with the default control,
// for the hat basis. Step halving is reproduced only where glm.fit() errors
// (no previous coefficients); a later halving is not reproduced.
Outcome glm_fit(int family, int link, const double* y, const double* offset,
                const double* prior, int n, const std::vector<int>& cell,
                const std::vector<double>& u, int m, const std::vector<double>& mustart,
                Fit& out, bool& conv) {
  auto predictor = [&](const std::vector<double>& coef, std::vector<double>& eta) {
    eta.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
      const int k = cell[i];
      const double lin = k < 0 ? 0.0 : (1.0 - u[i]) * coef[k] + u[i] * coef[k + 1];
      eta[i] = lin + offset[i];
    }
  };
  std::vector<double> eta(n), mu = mustart, z(n), w2(n), mev(n), coef;
  for (int i = 0; i < n; ++i) eta[i] = linkfun(link, mu[i]);
  // glm.fit() stops when the starting values are invalid.
  if (!valideta(link, eta) || !validmu(family, mu)) return ERROR;
  for (double e : eta) if (ISNAN(e)) return UNSUPPORTED;
  double devold = deviance(family, y, mu, prior, n), dev = devold;
  bool have_coef = false;
  conv = false;
  for (int iter = 1; iter <= 25; ++iter) {
    bool any_good = false;
    for (int i = 0; i < n; ++i) {
      mev[i] = mu_eta(link, eta[i]);
      if (!(prior[i] > 0.0)) { w2[i] = 0.0; z[i] = 0.0; continue; }
      // glm.fit() stops on NA or zero variances and NA derivatives.
      const double v = variance(family, mu[i]);
      if (ISNAN(v) || v == 0.0 || ISNAN(mev[i])) return ERROR;
      if (mev[i] == 0.0) { w2[i] = 0.0; z[i] = 0.0; continue; }
      z[i] = (eta[i] - offset[i]) + (y[i] - mu[i]) / mev[i];
      const double wi = std::sqrt((prior[i] * (mev[i] * mev[i])) / v);
      w2[i] = wi * wi;
      any_good = true;
    }
    if (!any_good) return UNSUPPORTED;
    if (!hat_wls(cell, u, z, w2, m, coef)) return UNSUPPORTED;
    predictor(coef, eta);
    for (int i = 0; i < n; ++i) mu[i] = linkinv(link, eta[i]);
    dev = deviance(family, y, mu, prior, n);
    if (!R_FINITE(dev) || !valideta(link, eta) || !validmu(family, mu)) {
      // glm.fit() halves the step towards the previous coefficients, or errors without them.
      return have_coef ? UNSUPPORTED : ERROR;
    }
    if (std::fabs(dev - devold) / (std::fabs(dev) + 0.1) < 1e-8) { conv = true; break; }
    devold = dev;
    have_coef = true;
  }
  out.theta = coef;
  out.deviance = dev;
  out.residuals.assign(n, 0.0);
  for (int i = 0; i < n; ++i) out.residuals[i] = (y[i] - mu[i]) / mu_eta(link, eta[i]);
  // glm.fit() returns the working weights of its last iteration's solve.
  out.weights = w2;
  return FIT;
}

// The family's default starting means, as family$initialize.
std::vector<double> cold_start(int family, const double* y, const double* prior, int n) {
  std::vector<double> mu(n);
  for (int i = 0; i < n; ++i) {
    switch (family) {
    case POISSON: mu[i] = y[i] + 0.1; break;
    case BINOMIAL: mu[i] = (prior[i] * (prior[i] == 0.0 ? 0.0 : y[i]) + 0.5) / (prior[i] + 1.0); break;
    default: mu[i] = y[i]; break;
    }
  }
  return mu;
}

// SplineReg_GLM(fast = TRUE, inits = guess): a warm start from the previous
// model's coefficients, and a cold restart from the family's default start if
// the warm fit fails, does not converge or ends above its starting deviance;
// the cold fit is kept only when its deviance is lower. Returns false when R
// would raise an error or take a branch not reproduced here.
bool irls_warm(int family, int link, const double* y, const double* offset,
               const double* prior, int n, const std::vector<int>& cell,
               const std::vector<double>& u, int m, const std::vector<double>& guess,
               Fit& out) {
  std::vector<double> mustart(n);
  for (int i = 0; i < n; ++i) {
    const int k = cell[i];
    const double lin = k < 0 ? 0.0 : (1.0 - u[i]) * guess[k] + u[i] * guess[k + 1];
    mustart[i] = linkinv(link, lin + offset[i]);
  }
  const double start_dev = deviance(family, y, mustart, prior, n);
  bool warm_conv = false;
  const Outcome warm = glm_fit(family, link, y, offset, prior, n, cell, u, m, mustart, out, warm_conv);
  if (warm == UNSUPPORTED) return false;
  const bool needs_cold = warm == ERROR || !warm_conv || !R_FINITE(out.deviance) ||
    out.deviance > start_dev + 1e-8 * std::max(1.0, std::fabs(start_dev));
  if (!needs_cold) return true;
  // family$initialize rejects responses outside the family's support; R then
  // keeps or rejects the warm fit through its own error handling.
  for (int i = 0; i < n; ++i) {
    const bool invalid = (family == POISSON && y[i] < 0.0) ||
      ((family == GAMMA || family == INVERSE_GAUSSIAN) && !(y[i] > 0.0)) ||
      (family == BINOMIAL && (y[i] < 0.0 || y[i] > 1.0)) ||
      (family == GAUSSIAN && link != IDENTITY && !(y[i] > 0.0));
    if (invalid) return false;
  }
  Fit cold;
  bool cold_conv = false;
  const Outcome cold_outcome = glm_fit(family, link, y, offset, prior, n, cell, u, m,
                                       cold_start(family, y, prior, n), cold, cold_conv);
  if (cold_outcome == UNSUPPORTED) return false;
  const bool prefer_cold = cold_outcome == FIT && R_FINITE(cold.deviance) &&
    (warm != FIT || !R_FINITE(out.deviance) || cold.deviance < out.deviance);
  if (prefer_cold) { out = cold; return true; }
  return warm == FIT;
}

}  // namespace

// Stage A of GenUnivariateFitter() in GD path mode (no parametric terms,
// finite plain boundaries, a standard family and link), continuing from the
// fit at iteration `j`: knot placement, coefficient guess for the new knot,
// warm-started IRLS fit, deviance, and the knot cap. Status: 0 knot cap
// reached, 1 iteration limit reached, 2 knot insertion failed, 3 the R code
// must refit iteration `j` from the returned knots and guess, 4 the R code must
// place the next knot after the returned fit at iteration `j`.
// R reference: the Stage A loop of GenUnivariateFitter() in R/UnivariateFitter.R.
// [[Rcpp::export]]
List gd_glm_path(NumericVector X, NumericVector Y, NumericVector offset,
                 NumericVector weights, LogicalVector positive,
                 IntegerVector site_counts, bool repeated_sites, NumericVector sites,
                 double beta, NumericVector extr, double tol, int family, int link,
                 NumericVector intknots_start, NumericVector guess_start,
                 NumericVector residuals_start, NumericVector working_start,
                 int j, double max_intknots, double last_iter) {
  const int n = X.size();
  const int m_sites = sites.size();
  std::vector<double> intknots(intknots_start.begin(), intknots_start.end());
  std::sort(intknots.begin(), intknots.end());
  std::vector<double> guess(guess_start.begin(), guess_start.end());
  std::vector<double> res(residuals_start.begin(), residuals_start.end());
  std::vector<double> ww(working_start.begin(), working_start.end());
  std::vector<double> loss, t;
  std::vector<int> cell;
  std::vector<double> u;
  int status;

  for (;;) {
    // Weighted residuals at the distinct positive-weight sites.
    NumericVector res_w(m_sites);
    if (repeated_sites) {
      int i = 0;
      for (int g = 0; g < m_sites; ++g) {
        long double rsum = 0.0, wsum = 0.0;
        for (int c = 0; c < site_counts[g]; ++i) {
          if (positive[i] != TRUE) continue;
          rsum += res[i];
          wsum += weights[i] * ww[i];
          ++c;
        }
        res_w[g] = ((double) rsum / (double) site_counts[g]) * (double) wsum;
      }
    } else {
      for (int i = 0, g = 0; i < n; ++i) {
        if (positive[i] == TRUE) res_w[g++] = res[i] * weights[i] * ww[i];
      }
    }
    bool finite = true;
    for (int g = 0; g < m_sites; ++g) if (!R_FINITE(res_w[g])) finite = false;
    // Same-sign clusters, as rle(sign(res.weighted)).
    std::vector<int> ends;
    for (int g = 1; g < m_sites; ++g) {
      const double a = res_w[g - 1], b = res_w[g];
      if (((a > 0) - (a < 0)) != ((b > 0) - (b < 0))) ends.push_back(g);
    }
    ends.push_back(m_sites);
    const int nc = ends.size();
    std::vector<double> cs(m_sites);
    long double running = 0.0;
    for (int g = 0; g < m_sites; ++g) { running += res_w[g]; cs[g] = (double) running; }
    NumericVector dcum(nc), cw(nc);
    std::vector<double> means(nc), ranges(nc);
    double max_mean = 0.0, max_range = 0.0;
    for (int c = 0; c < nc; ++c) {
      const int end = ends[c], begin = c == 0 ? 0 : ends[c - 1];
      dcum[c] = end;
      const double previous = c == 0 ? 0.0 : cs[begin - 1];
      means[c] = std::fabs((cs[end - 1] - previous) / (double) (end - begin));
      ranges[c] = sites[end - 1] - sites[begin];
      max_mean = std::max(max_mean, means[c]);
      max_range = std::max(max_range, ranges[c]);
    }
    if (!finite || !(max_mean > 0.0) || !R_FINITE(max_mean)) { status = 4; break; }
    for (int c = 0; c < nc; ++c) {
      const double range = max_range != 0.0 ? ranges[c] / max_range : ranges[c];
      cw[c] = beta * (means[c] / max_mean) + (1.0 - beta) * range;
    }
    NumericVector oldknots(4 + intknots.size());
    oldknots[0] = extr[0]; oldknots[1] = extr[1]; oldknots[2] = extr[0]; oldknots[3] = extr[1];
    for (size_t k = 0; k < intknots.size(); ++k) oldknots[4 + k] = intknots[k];
    const double newknot = Knotnew(cw, res_w, sites, dcum, oldknots, tol, 2)[0];
    if (!R_FINITE(newknot) || newknot <= extr[0] || newknot >= extr[1]) { status = 2; break; }

    // newknot.guess(): the current spline's value at the new knot, inserted
    // after the coefficients of the knots below it.
    t.assign(1, extr[0]);
    t.insert(t.end(), intknots.begin(), intknots.end());
    t.push_back(extr[1]);
    std::vector<int> nk_cell; std::vector<double> nk_u;
    hat_basis(&newknot, 1, t, nk_cell, nk_u);
    const int position = std::count_if(intknots.begin(), intknots.end(),
                                       [&](double k) { return k < newknot; });
    double value = 0.0;
    if (nk_cell[0] >= 0) {
      const int k = nk_cell[0];
      long double s = (long double) ((1.0 - nk_u[0]) * guess[k]) + (long double) (nk_u[0] * guess[k + 1]);
      value = (double) s;
    }
    guess.insert(guess.begin() + position + 1, value);
    intknots.insert(intknots.begin() + position, newknot);

    if ((double) j >= last_iter) { status = 1; break; }
    ++j;

    t.assign(1, extr[0]);
    t.insert(t.end(), intknots.begin(), intknots.end());
    t.push_back(extr[1]);
    for (size_t k = 1; k < t.size(); ++k) if (!(t[k] > t[k - 1])) { status = 3; goto done; }
    hat_basis(X.begin(), n, t, cell, u);
    {
      Fit fit;
      if (!irls_warm(family, link, Y.begin(), offset.begin(), weights.begin(), n, cell, u,
                     t.size(), guess, fit)) {
        status = 3; break;
      }
      guess = fit.theta;
      res = fit.residuals;
      ww = fit.weights;
      loss.push_back(fit.deviance);
    }
    if ((double) j == max_intknots + 1.0) { status = 0; break; }
  }
done:
  return List::create(Named("status") = status, Named("j") = j,
                      Named("intknots") = NumericVector(intknots.begin(), intknots.end()),
                      Named("guess") = NumericVector(guess.begin(), guess.end()),
                      Named("loss") = NumericVector(loss.begin(), loss.end()),
                      Named("residuals") = NumericVector(res.begin(), res.end()),
                      Named("weights") = NumericVector(ww.begin(), ww.end()));
}
