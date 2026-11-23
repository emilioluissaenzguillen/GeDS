/* Routines for GeDS R-package
 *
 * Routines are based on the algorithm in Kaishev et al. (2016).
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 *
 * These functions are distributed WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.	 See the GNU General Public License for more details.
 */


#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
int whmx(NumericVector vector) {
  int i, k=0, length = vector.size();
  double max = vector(0);
  for(i=0; i<length; i++){
    if(vector(i) > max) {
      max = vector(i);
      k = i;
    }
  }
  return k;
}

#include <Rcpp.h>
#include <algorithm>

// [[Rcpp::export]]
bool is_internal_knot(double newknot,
                             Rcpp::NumericVector sortedknots,
                             double tol) {
  if (sortedknots.size() == 0) return false;
  
  double min_k = *std::min_element(sortedknots.begin(), sortedknots.end());
  double max_k = *std::max_element(sortedknots.begin(), sortedknots.end());
  
  bool newknot_is_internal =
    (newknot >= min_k + tol) &&
    (newknot <= max_k - tol);
  
  return newknot_is_internal;
}

namespace {
bool has_ordered_sw_support(const NumericVector& knots,
                            const NumericVector& sites, int order, double tol) {
  int next_site = 0;
  int n_basis = knots.size() - order;
  double previous_site = R_NegInf;
  for (int i = 0; i < n_basis; ++i) {
    bool matched = false;
    while (next_site < sites.size()) {
      double z = sites[next_site++];
      if (z <= previous_site) continue;
      bool interior = knots[i] + tol < z && z < knots[i + order] - tol;
      bool left_endpoint = i == 0 && z == knots[0] && knots[order - 1] == knots[0];
      bool right_endpoint = i == n_basis - 1 && z == knots[knots.size() - 1] &&
        knots[n_basis] == knots[knots.size() - 1];
      if (interior || left_endpoint || right_endpoint) {
        previous_site = z;
        matched = true;
        break;
      }
    }
    if (!matched) return false;
  }
  return true;
}

// Same greedy matching for strictly increasing sites. Each basis takes the
// first unused site satisfying one of its conditions, so the scan can jump to
// the first candidate of each condition by binary search.
bool has_ordered_sw_support_sorted(const NumericVector& knots,
                                   const NumericVector& sites, int order, double tol) {
  const double *s = sites.begin(), *end = sites.end();
  const int m = sites.size(), n_knots = knots.size(), n_basis = n_knots - order;
  int next_site = 0;
  for (int i = 0; i < n_basis; ++i) {
    int match = m;
    const double lower = knots[i] + tol, upper = knots[i + order] - tol;
    int j = (int)(std::upper_bound(s + next_site, end, lower) - s);
    if (j < m && s[j] < upper) match = j;
    if (i == 0 && knots[order - 1] == knots[0]) {
      j = (int)(std::lower_bound(s + next_site, end, (double)knots[0]) - s);
      if (j < m && s[j] == knots[0] && j < match) match = j;
    }
    if (i == n_basis - 1 && knots[n_basis] == knots[n_knots - 1]) {
      j = (int)(std::lower_bound(s + next_site, end, (double)knots[n_knots - 1]) - s);
      if (j < m && s[j] == knots[n_knots - 1] && j < match) match = j;
    }
    if (match == m) return false;
    next_site = match + 1;
  }
  return true;
}

bool strictly_increasing(const NumericVector& v) {
  for (int i = 1; i < v.size(); ++i) {
    if (!(v[i] > v[i - 1])) return false;
  }
  return true;
}
}

// R reference: Knotnew_R_clean() (Knotnew_R() is the original) in R/cpp_functions_R.R.
// [[Rcpp::export]]
NumericVector Knotnew(NumericVector weights, NumericVector residuals, NumericVector x,
                      NumericVector dcum, NumericVector oldknots, double tol,
                      int support_order = 2) {
  
  int u = dcum.size();
  int n_oldknots = oldknots.size();
  int n_oldintknots = n_oldknots - 2 * support_order;
  int data_size = x.size();

  if (u == 0 || weights.size() == 0) {
    return NumericVector::create(NA_REAL, NA_REAL);
  }

  // Callers pass distinct sorted sites; other inputs keep the linear scans.
  const bool sorted_sites = strictly_increasing(x);
  NumericVector sw_sites = x;
  if (!sorted_sites) {
    sw_sites = clone(x);
    std::sort(sw_sites.begin(), sw_sites.end());
  }
  std::vector<double> oldintknots(n_oldintknots > 0 ? n_oldintknots : 0);
  for (int k = 0; k < n_oldintknots; ++k) {
    oldintknots[k] = oldknots[2 * support_order + k];
  }

  NumericVector candidate_weights = clone(weights);
  
  int best_index = 0, dcumInf, dcumSup, i, j;
  long double sup, inf;
  long double newknot = NA_REAL;
  bool has_old_knots_between, is_valid_knot;
  bool knot_found = false;
  
  // Iterate through each cluster to find the best knot position
  for (int cluster = 0; cluster < u; ++cluster) {
    
    // Find the index of the cluster with the highest weight
    best_index = whmx(candidate_weights);
    
    // If all weights are zero, fall back to triangular profile and recompute best_index
    if (candidate_weights[best_index] == 0.0) {
      int n = static_cast<int>(candidate_weights.size());
      double center = (n + 1.0) / 2.0;
      // Reassign triangular weights
      for (int i = 0; i < n; ++i) {
        candidate_weights[i] = 1.0 - std::abs((double(i + 1) - center)) / center;
      }
      // Recompute best_index (mimic R's which.max: first occurrence)
      best_index = std::distance(
        candidate_weights.begin(),
        std::max_element(candidate_weights.begin(), candidate_weights.end())
      );
    }
    
    // Determine lower and upper bounds
    dcumInf = (best_index == 0) ? 0 : dcum[best_index - 1];
    dcumSup = dcum[best_index] - 1;
    
    // Compute superior and inferior x bounds
    sup = x(dcumSup);
    inf = x(dcumInf);
    
    // 1) If there are previous internal knots, check whether they lie within the current cluster
    has_old_knots_between = false;
    
    if (n_oldintknots > 0) {
      if (inf != sup) {
        // For each internal knot, check if it lies within the current cluster's boundaries
        for (j = 0; j < n_oldintknots; ++j) {
          has_old_knots_between = has_old_knots_between || 
            ((inf <= oldintknots[j]) && (sup >= oldintknots[j]));
        }
      } else { 
        // If inf == sup (single point), check if the singleton is close to previous knots
        for (j = 0; j < n_oldintknots; ++j) {
          has_old_knots_between = has_old_knots_between || 
            (std::abs(inf - oldintknots[j]) < tol);
        }  
      }  
    } 
    
    // 2) If an internal knot already exists in this cluster, set its weight to zero
    if (has_old_knots_between) {
      candidate_weights[best_index] = 0;
      continue;
    }     
    
    // 3) Compute the new knot as a weighted average of x-coordinates in the current cluster
    newknot = std::inner_product(
      residuals.begin() + dcumInf, residuals.begin() + dcumSup + 1,
      x.begin() + dcumInf, 0.0) /
        std::accumulate(residuals.begin() + dcumInf, residuals.begin() + dcumSup + 1, 0.0);   
    
    // Insert the new knot and sort all knots
    NumericVector sortedknots(n_oldknots + 1);
    sortedknots[n_oldknots] = newknot;
    
    for (i = 0; i < n_oldknots; ++i) {
      sortedknots[i] = oldknots[i];
    }     
    
    std::sort(sortedknots.begin(), sortedknots.end());
    
    // Step 4: Retain the minimum-support rule and require distinct-site SW matching.
    is_valid_knot = true;
    
    for (i = 0; i < n_oldknots - (support_order - 1); ++i) {
      bool valid_interval = false;

      if (sorted_sites) {
        // The first site above the left bound is the only one to check.
        const double lower = sortedknots[i] + tol;
        j = (int)(std::upper_bound(x.begin(), x.end(), lower) - x.begin());
        valid_interval = j < data_size && x[j] < sortedknots[i + support_order] - tol;
      } else {
        for (j = 0; j < data_size; ++j) {
          valid_interval = valid_interval ||
            ((sortedknots[i] + tol < x[j]) &&
             (x[j] < sortedknots[i + support_order] - tol));
          if (valid_interval) break;
        }
      }

      is_valid_knot = is_valid_knot && valid_interval;
      if (!is_valid_knot) break;
    }

    if (is_valid_knot) {
      is_valid_knot = sorted_sites ?
        has_ordered_sw_support_sorted(sortedknots, sw_sites, support_order, tol) :
        has_ordered_sw_support(sortedknots, sw_sites, support_order, tol);
    }
    
    // 5) Check new knot is not a boundary knot
    bool newknot_is_internal = is_internal_knot(newknot, sortedknots, tol);
    
    // If the newknot placement is invalid, set the weight of the cluster to zero
    if (!is_valid_knot || !newknot_is_internal) {
      candidate_weights[best_index] = 0;
    } else { 
      knot_found = true;
      break; // Break out of the loop if a valid knot is found
    }  
  }   

  if (!knot_found) {
    return NumericVector::create(NA_REAL, NA_REAL);
  }
  
  // Return the new knot and its index (adjusted to 1-based for R)
  return NumericVector::create(newknot, best_index + 1);
  
}      


// [[Rcpp::export]]
NumericVector makenewknots(NumericVector knots, int degree) {
  int i, j, length;
  double cumtemp;
  length = knots.size();
  NumericVector newknots(length-(degree-2));
  for(i=0; i< (length-(degree-2)); ++i){
    for(j=0, cumtemp = 0.0; j<(degree-1); ++j) cumtemp += knots[i+j];
    newknots(i) = cumtemp / (degree - 1);
  }
  return newknots;
}

// [[Rcpp::export]]
NumericVector makeEpsilonsb(NumericVector data, NumericVector Xs, NumericVector Ys, int degree) {
  int length = Xs.size();
  int i, j, p = data.size();
  double cumtemp;
  p = length - degree;
  NumericVector epsilon(p);
  for(i=0; i<p ; ++i){
    for(j=1, cumtemp = 0.0; j<degree; ++j) cumtemp += Xs[i+j];
    epsilon(i) = cumtemp / (degree - 1);
    }
  return epsilon;
  }

NumericVector ctrlpolyfun(NumericVector data, NumericVector Xs, NumericVector Ys, int degree) {
  int length = Xs.size();
  int i, j, p = length - degree, n = data.size();

  if (p < 2) {
    stop("ctrlpolyfun requires at least two epsilon values");
  }

  NumericVector epsilons = makeEpsilonsb(data, Xs, Ys, degree);
  NumericVector vector(n);
    for (i=0; i<n; i++){
      // Use the final interval for values beyond the last epsilon.
      int k = p - 1;
      for (j = 1; j < p; j++) {
        if (epsilons[j] >= data[i]){
          k =  j; break;
          }
        }
      vector[i] =Ys(k - 1) + (Ys(k) - Ys(k - 1)) * (data(i) - epsilons(k - 1)) / (epsilons(k) - epsilons(k - 1));
      }
    return vector;
    }

// R reference: tensorProd_R() in R/cpp_functions_R.R.
// [[Rcpp::export]]
NumericMatrix tensorProd(NumericMatrix Xmat, NumericMatrix Ymat){
  int ii, jjx, jjy, kk, ib = Xmat.nrow(), jx = Xmat.ncol(), jy = Ymat.ncol();
  NumericMatrix ris(ib,jx*jy);
  for(jjx=0, kk=0; jjx<jx; jjx++){
    for(jjy=0;jjy<jy; jjy++){
      for(ii=0; ii<ib; ii++){
        ris(ii,kk) = Xmat(ii,jjx)*Ymat(ii,jjy);
      }
      kk++;
    }
  }
  return ris;
}

// [[Rcpp::export]]
NumericMatrix makeNewMatrCPP(NumericMatrix matrix, Nullable<IntegerMatrix> tab = R_NilValue, bool by_row = false) {
  int num_rows = matrix.nrow();
  
  // Handle NULL `tab`
  if (tab.isNull()) {
    return matrix;
  } 
  
  IntegerMatrix tab_mat(tab);
  int tab_rows = tab_mat.nrow();
  int tab_cols = tab_mat.ncol();
  
  std::vector<int> recurr;
  
  // Flatten `tab`
  if (by_row) {
    for (int i = 0; i < tab_rows; ++i) {
      for (int j = 0; j < tab_cols; ++j) {
        recurr.push_back(tab_mat(i, j));  // Mimics `c(t(tab))`
      } 
    }
  } else {
    for (int j = 0; j < tab_cols; ++j) {
      for (int i = 0; i < tab_rows; ++i) {
        recurr.push_back(tab_mat(i, j));  // Mimics `c(tab)`
      }
    } 
  }
  
  // Remove zeros
  recurr.erase(std::remove(recurr.begin(), recurr.end(), 0), recurr.end());
   
  // Compute cumulative sum
  std::vector<int> ids(recurr.size() + 1, 0);
  for (size_t i = 0; i < recurr.size(); ++i) {
    ids[i + 1] = ids[i] + recurr[i];
  }
   
  int n = ids.size() - 1;
  NumericMatrix ret(n, 3);
   
  // Compute newX, newY, and newres
  for (int i = 0; i < n; ++i) {
    if (ids[i] >= num_rows) break; // Prevent out-of-bounds access
     
    double sum_res = 0;
    for (int j = ids[i]; j < std::min(ids[i+1], num_rows); ++j) {
      sum_res += matrix(j, 2);
    }
     
    // FIX: Take X and Y from row `i`, not `ids[i]`
    ret(i, 0) = matrix(i, 0); // FIXED: Now `newX`
    ret(i, 1) = matrix(i, 1); // FIXED: Now `newY`
    ret(i, 2) = sum_res;      // Sum correctly assigned
  } 
  
  return ret;
} 

// R reference: findNewDimKnot_R() in R/cpp_functions_R.R.
// [[Rcpp::export]]
List findNewDimKnot(
    IntegerVector dcumFixedDim_Dim,
    NumericVector Dim_weights,
    NumericVector Dim_oldknots,
    NumericVector Dim_values,
    NumericVector residuals)
{
  int u = dcumFixedDim_Dim.size(); // Number of clusters
  bool flagDim = false; // Flag for negative weights
  double Dim_newknot = NumericVector::get_na();  
  double weightDim = NumericVector::get_na();    
  int dcumInf = -1, dcumSup = -1;  // Initialize with invalid values

  if (Dim_weights.size() != u) {
    stop("Dim_weights and dcumFixedDim_Dim must have the same length.");
  }
  if (Dim_values.size() != residuals.size()) {
    stop("Dim_values and residuals must have the same length.");
  }

  // Work on a local copy so rejected clusters do not mutate the R input.
  NumericVector candidate_weights = clone(Dim_weights);
  
  for (int i = 0; i < u; ++i) {
    if (is_true(all(candidate_weights < 0))) {
      flagDim = true;
      break;
    }
    
    // Find the index of the cluster with the highest weight
    int best_index = which_max(candidate_weights);
    
    if (best_index < 0 || best_index >= dcumFixedDim_Dim.size()) {
      stop("Invalid cluster index detected.");
    } 
    
    // Determine boundaries indexes
    dcumInf = (best_index == 0) ? 0 : dcumFixedDim_Dim[best_index - 1];
    dcumSup = dcumFixedDim_Dim[best_index] - 1;
    
    if (dcumInf >= Dim_values.size() || dcumSup >= Dim_values.size() || dcumInf > dcumSup) {
      stop("Indexing error: dcumInf or dcumSup out of range.");
    } 
    
    // Calculate superior and inferior Dim-bounds
    double sup = Dim_values[dcumSup];
    double inf = Dim_values[dcumInf];
    
    // Compute the residual-weighted average of the coordinate values.
    Dim_newknot = std::inner_product(
      residuals.begin() + dcumInf, residuals.begin() + dcumSup + 1,
      Dim_values.begin() + dcumInf,
      0.0
    ) / std::accumulate(
      residuals.begin() + dcumInf,
      residuals.begin() + dcumSup + 1,
      0.0
    );
    
    // Validate new knot conditions
    bool cond1 = (dcumSup - dcumInf) != 0;
    bool cond2 = !is_true(any((Dim_oldknots >= inf) & (Dim_oldknots <= sup)));
    bool cond3 = !is_true(any(abs(inf - Dim_oldknots) < 1e-12));
    
    if ((cond1 && cond2) || (!cond1 && cond3)) {
      weightDim = candidate_weights[best_index];
      break;
    } else {  
      candidate_weights[best_index] = R_NegInf;
    }
  }
  
  // Check if all Dim.weights were turned to -Inf
  if (is_true(all(candidate_weights < 0))) {
    flagDim = true;
    }
  
  return List::create(
    Named("Dim.newknot") = Dim_newknot,
    Named("weightDim") = weightDim,
    Named("flagDim") = flagDim,
    Named("dcumInf") = dcumInf + 1,  // Convert to 1-based indexing
    Named("dcumSup") = dcumSup + 1  // Convert to 1-based indexing
  );
} 




// Weighted least-squares fit of a linear spline with breakpoints t (boundary
// knots first and last, strictly increasing). The hat-function basis makes the
// normal equations tridiagonal, so they are solved by an LDL' sweep in O(n).
// Matches splineDesign(ord = 2, outer.ok = TRUE): points outside [t0, tm] have
// a zero basis row. Returns NULL when the system is singular or poorly
// conditioned so the caller can fall back to a pivoted QR.
namespace {
bool linear_spline_ls_core(const double* x, const double* y, const double* w, int n,
                           const double* t, int m, std::vector<double>& theta,
                           std::vector<double>& fitted) {
  if (m < 2) return false;
  for (int k = 1; k < m; k++) {
    if (!(t[k] > t[k - 1])) return false;
  }
  std::vector<double> d(m, 0.0), e(m - 1, 0.0), b(m, 0.0), u(n, 0.0);
  std::vector<int> cell(n, -1);
  const double *tb = t, *te = t + m;
  for (int i = 0; i < n; i++) {
    const double xi = x[i];
    if (xi < t[0] || xi > t[m - 1]) continue;
    // Interval k with t[k] <= x < t[k + 1]; the right boundary closes the last.
    int k = (int)(std::upper_bound(tb, te, xi) - tb) - 1;
    if (k > m - 2) k = m - 2;
    const double ui = (xi - t[k]) / (t[k + 1] - t[k]), ai = 1.0 - ui;
    cell[i] = k;
    u[i] = ui;
    const double wi = w[i];
    if (wi == 0.0) continue;
    d[k] += wi * ai * ai;
    d[k + 1] += wi * ui * ui;
    e[k] += wi * ai * ui;
    b[k] += wi * ai * y[i];
    b[k + 1] += wi * ui * y[i];
  }
  const double rel_tol = 1e-10;
  if (!(d[0] > 0.0)) return false;
  double pivot = d[0];
  std::vector<double> l(m, 0.0);
  for (int k = 1; k < m; k++) {
    l[k] = e[k - 1] / pivot;
    const double next = d[k] - l[k] * e[k - 1];
    if (!(d[k] > 0.0) || !(next > rel_tol * d[k])) return false;
    b[k] -= l[k] * b[k - 1];
    pivot = d[k] = next;
  }
  theta.assign(m, 0.0);
  theta[m - 1] = b[m - 1] / d[m - 1];
  for (int k = m - 2; k >= 0; k--) {
    theta[k] = (b[k] - e[k] * theta[k + 1]) / d[k];
  }
  fitted.assign(n, 0.0);
  for (int i = 0; i < n; i++) {
    const int k = cell[i];
    if (k >= 0) fitted[i] = (1.0 - u[i]) * theta[k] + u[i] * theta[k + 1];
  }
  return true;
}
}

// R reference: .linear_spline_ls_R() in R/cpp_functions_R.R.
// [[Rcpp::export]]
SEXP linear_spline_ls(NumericVector x, NumericVector y, NumericVector w,
                      NumericVector t) {
  const int n = x.size();
  if (y.size() != n || w.size() != n) return R_NilValue;
  std::vector<double> theta, fitted;
  if (!linear_spline_ls_core(x.begin(), y.begin(), w.begin(), n, t.begin(), t.size(),
                             theta, fitted)) {
    return R_NilValue;
  }
  return List::create(Named("theta") = NumericVector(theta.begin(), theta.end()),
                      Named("fitted") = NumericVector(fitted.begin(), fitted.end()));
}

// One Stage A knot-placement step of the univariate Gaussian fitter
// (Kaishev et al., 2016, steps 2-7). It reproduces the R implementation
// operation by operation: weighted residuals at positive-weight sites,
// averaged over repeated sites with R's long-double sum(), numerical zeros
// removed, same-sign clusters, normalised cluster means and ranges, and
// Knotnew(). Returns NA when no valid knot exists or all residuals vanish.
// R reference: .geds_knot_step_R() in R/cpp_functions_R.R.
// [[Rcpp::export]]
double geds_knot_step(NumericVector residuals, LogicalVector positive,
                      NumericVector placement_weights, IntegerVector site_counts,
                      bool repeated_sites, NumericVector site_scale,
                      double geometry_tol, NumericVector sites, double beta,
                      NumericVector extr, NumericVector intknots, double tol,
                      int support_order) {
  const int n_pos = placement_weights.size();
  std::vector<double> contribution(n_pos);
  for (int i = 0, k = 0; i < residuals.size(); ++i) {
    if (positive[i] != TRUE) continue;
    if (k >= n_pos) stop("Inconsistent positive-weight sites.");
    contribution[k] = residuals[i] * placement_weights[k];
    ++k;
  }

  const int m = sites.size();
  NumericVector res_w(m);
  if (repeated_sites) {
    if (site_counts.size() != m) stop("Inconsistent repeated-site counts.");
    for (int g = 0, start = 0; g < m; ++g) {
      long double total = 0.0;
      for (int k = start; k < start + site_counts[g]; ++k) total += contribution[k];
      res_w[g] = (double) total / (double) site_counts[g];
      start += site_counts[g];
    }
  } else {
    if (n_pos != m) stop("Inconsistent placement sites.");
    for (int g = 0; g < m; ++g) res_w[g] = contribution[g];
  }

  bool all_zero = true;
  for (int g = 0; g < m; ++g) {
    if (std::fabs(res_w[g]) <= geometry_tol * site_scale[g]) res_w[g] = 0.0;
    if (!(res_w[g] == 0.0)) all_zero = false;
  }
  if (all_zero || m == 0) return NA_REAL;

  // Same-sign clusters, as rle(sign(res_w)); NaN never joins a cluster.
  std::vector<double> sign(m);
  for (int g = 0; g < m; ++g) {
    const double v = res_w[g];
    sign[g] = std::isnan(v) ? v : (v > 0.0) - (v < 0.0);
  }
  std::vector<int> ends;
  for (int g = 1; g < m; ++g) {
    if (!(sign[g] == sign[g - 1])) ends.push_back(g);
  }
  ends.push_back(m);
  const int u = ends.size();

  // Cluster means from R's cumsum(), which accumulates in long double.
  std::vector<double> cs(m);
  long double running = 0.0;
  for (int g = 0; g < m; ++g) {
    running += res_w[g];
    cs[g] = (double) running;
  }
  NumericVector dcum(u), weights(u);
  std::vector<double> means(u), ranges(u);
  double max_mean = R_NegInf, max_range = R_NegInf;
  bool mean_nan = false, range_nan = false;
  for (int c = 0; c < u; ++c) {
    const int end = ends[c], begin = c == 0 ? 0 : ends[c - 1];
    dcum[c] = end;
    const double previous = c == 0 ? 0.0 : cs[begin - 1];
    means[c] = std::fabs((cs[end - 1] - previous) / (double) (end - begin));
    ranges[c] = sites[end - 1] - sites[begin];
    if (std::isnan(means[c])) mean_nan = true; else if (means[c] > max_mean) max_mean = means[c];
    if (std::isnan(ranges[c])) range_nan = true; else if (ranges[c] > max_range) max_range = ranges[c];
  }
  if (mean_nan) max_mean = NA_REAL;
  if (range_nan) stop("Missing cluster ranges in knot placement.");
  for (int c = 0; c < u; ++c) {
    const double mean = means[c] / max_mean;
    const double range = max_range != 0.0 ? ranges[c] / max_range : ranges[c];
    weights[c] = beta * mean + (1.0 - beta) * range;
  }

  const int n_ext = extr.size();
  NumericVector oldknots(n_ext * support_order + intknots.size());
  for (int r = 0, k = 0; r < support_order; ++r) {
    for (int e = 0; e < n_ext; ++e) oldknots[k++] = extr[e];
  }
  for (int k = 0; k < intknots.size(); ++k) oldknots[n_ext * support_order + k] = intknots[k];

  return Knotnew(weights, res_w, sites, dcum, oldknots, tol, support_order)[0];
}

// Empirical semivariance of residuals over nearest neighbours, as
// colMeans((r - matrix(r[nn], nrow = length(r), byrow = TRUE))^2) / 2, where
// nn is the p x n matrix of neighbour indices (1-based). Matches R's colMeans,
// which accumulates in long double and divides before rounding.
// R reference: .gd_neighbour_semivariance_R() in R/cpp_functions_R.R.
// [[Rcpp::export]]
NumericVector gd_neighbour_semivariance(NumericVector r, IntegerMatrix nn) {
  const int p = nn.nrow(), n = nn.ncol();
  if (n != r.size()) stop("'nn' must have one column per residual.");
  NumericVector out(p);
  for (int j = 0; j < p; ++j) {
    long double total = 0.0;
    for (int i = 0; i < n; ++i) {
      const double d = r[i] - r[nn(j, i) - 1];
      total += d * d;
    }
    total /= n;
    out[j] = (double) total / 2.0;
  }
  return out;
}

// Stage A of the univariate Gaussian fitter in GD path mode (no parametric
// terms, finite plain boundaries). Mirrors the R loop in UnivariateFitter():
// linear least-squares fit, weighted RSS as R's long-double sum(), stop at
// the knot cap, knot placement (geds_knot_step) and the boundary check.
// `first_residuals` replaces the first fit when a starting fit or a staged
// restart supplies it; `append_first_rss` adds its RSS (starting fit only).
// Returns NULL when the tridiagonal solver declines a fit; the R loop then
// recomputes the whole path with its pivoted-QR fallback.
// R reference: the Stage A loop of UnivariateFitter() in R/UnivariateFitter.R.
// [[Rcpp::export]]
SEXP gd_stage_a_path(NumericVector X, NumericVector Y, NumericVector offset,
                     NumericVector weights, LogicalVector positive,
                     NumericVector placement_weights, IntegerVector site_counts,
                     bool repeated_sites, NumericVector site_scale, double geometry_tol,
                     NumericVector sites, double beta, NumericVector extr, double tol,
                     double max_intknots, NumericVector intknots_start,
                     NumericVector loss_start, Nullable<NumericVector> first_residuals,
                     bool append_first_rss, int first_iter, double last_iter) {
  const int n = X.size();
  if (Y.size() != n || offset.size() != n || weights.size() != n || extr.size() != 2) {
    return R_NilValue;
  }
  std::vector<double> y0(n);
  for (int i = 0; i < n; ++i) y0[i] = Y[i] - offset[i];
  std::vector<double> intknots(intknots_start.begin(), intknots_start.end());
  std::vector<double> loss(loss_start.begin(), loss_start.end());
  NumericVector residuals(n);
  std::vector<double> theta, fitted, full, inner;
  const double boundary_tolerance = std::sqrt(DBL_EPSILON);

  auto weighted_rss = [&](const NumericVector& r) {
    long double total = 0.0;
    for (int i = 0; i < n; ++i) total += weights[i] * (r[i] * r[i]);
    return (double) total;
  };

  int j = first_iter;
  for (; (double) j <= last_iter; ++j) {
    if (j == first_iter && first_residuals.isNotNull()) {
      NumericVector start(first_residuals.get());
      if (start.size() != n) return R_NilValue;
      for (int i = 0; i < n; ++i) residuals[i] = start[i];
      if (append_first_rss) loss.push_back(weighted_rss(residuals));
    } else {
      // Breakpoints as sort(c(intknots, rep(extr, 2))) without its two ends.
      full.assign(intknots.begin(), intknots.end());
      full.push_back(extr[0]); full.push_back(extr[1]);
      full.push_back(extr[0]); full.push_back(extr[1]);
      std::sort(full.begin(), full.end());
      inner.assign(full.begin() + 1, full.end() - 1);
      if (!linear_spline_ls_core(X.begin(), y0.data(), weights.begin(), n,
                                 inner.data(), inner.size(), theta, fitted)) {
        return R_NilValue;
      }
      for (int i = 0; i < n; ++i) residuals[i] = Y[i] - (fitted[i] + offset[i]);
      loss.push_back(weighted_rss(residuals));
    }
    if ((double) j == max_intknots + 1.0) break;

    NumericVector current(intknots.begin(), intknots.end());
    const double newknot = geds_knot_step(residuals, positive, placement_weights,
                                          site_counts, repeated_sites, site_scale,
                                          geometry_tol, sites, beta, extr, current,
                                          tol, 2);
    if (ISNAN(newknot)) break;
    if (R_FINITE(newknot)) {
      const double scale = std::fabs(newknot) > boundary_tolerance ? std::fabs(newknot) : 1.0;
      if (std::fabs(newknot - extr[0]) / scale <= boundary_tolerance ||
          std::fabs(newknot - extr[1]) / scale <= boundary_tolerance) break;
    }
    intknots.push_back(newknot);
  }
  // The R loop leaves j at its last value when it runs to completion.
  if ((double) j > last_iter) --j;

  return List::create(Named("loss") = NumericVector(loss.begin(), loss.end()),
                      Named("intknots") = NumericVector(intknots.begin(), intknots.end()),
                      Named("residuals") = residuals, Named("j") = j);
}
