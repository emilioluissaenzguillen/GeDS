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

// [[Rcpp::export]]
NumericVector Knotnew(NumericVector weights, NumericVector residuals, NumericVector x,
                      NumericVector dcum, NumericVector oldknots, double tol) {
  
  int u = dcum.size();
  int n_oldknots = oldknots.size();
  int n_oldintknots = n_oldknots - 6;
  int data_size = x.size();
  
  int best_index, dcumInf, dcumSup, i, j;
  long double sup, inf, newknot;
  bool has_old_knots_between, is_valid_knot;
  
  // Iterate through each cluster to find the best knot position
  for (int cluster = 0; cluster < u; ++cluster) {
    
    // Find the index of the cluster with the highest weight
    best_index = whmx(weights);
    
    // Determine lower and upper bounds
    dcumInf = (best_index == 0) ? 0 : dcum[best_index - 1];
    dcumSup = dcum[best_index] - 1;
    
    // Compute superior and inferior x bounds
    sup = x(dcumSup);
    inf = x(dcumInf);
    
    // 1) If there are previous internal knots, check whether they lie within the current cluster
    has_old_knots_between = false;
    
    if (n_oldintknots > 0) {
      NumericVector oldintknots(n_oldintknots);
      
      // Extract internal knots (ignoring the boundary knots)
      for (i = 0; i < n_oldintknots; ++i) {
        oldintknots[i] = oldknots[6 + i];
      }
      
      if (inf != sup) {
        // For each internal knot, check if it lies within the current cluster's boundaries
        for (j = 0; j < n_oldintknots; ++j) {
          has_old_knots_between = has_old_knots_between || 
            ((inf <= oldintknots[j]) && (sup >= oldintknots[j]));
        }
      } else { 
        // If inf == sup (single point), check if the knot is close to the boundary
        for (j = 0; j < n_oldintknots; ++j) {
          has_old_knots_between = has_old_knots_between || 
            (std::abs(inf - oldintknots[j]) < tol);
        }  
      }  
    } 
    
    // 2) If an internal knot already exists in this cluster, set its weight to zero
    if (has_old_knots_between) {
      weights[best_index] = 0;
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
    
    // Step 4: Check if new knot placement satisfies the minimum support constraint, i.e.,
    // For each consecutive set of 4 sortedknots, check whether there is at least one x
    // that falls between the 1st and 4th knot in that set
    is_valid_knot = true;
    
    for (i = 0; i < n_oldknots - 2; ++i) {
      bool valid_interval = false;
      
      for (j = 0; j < data_size; ++j) {
        valid_interval = valid_interval || 
          ((sortedknots[i] < x[j]) && (sortedknots[i + 3] > x[j])); 
        
        if (valid_interval) break;
      }     
      
      is_valid_knot = is_valid_knot && valid_interval;
      if (!is_valid_knot) break;
    }     
    
    // If the newknot placement is invalid, set the weight of the cluster to zero
    if (!is_valid_knot) {
      weights[best_index] = 0;
    } else { 
      // Break out of the loop if a valid knot is found
      break;
    }  
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
  NumericVector epsilons;
  int length = Xs.size();
  int i, j, k, p, n = data.size();
  epsilons = makeEpsilonsb(data, Xs, Ys, degree);
  NumericVector vector(n);
  p = length - degree;
    for (i=0; i<n; i++){
      for (j = 1; j < p; j++) {
        if (epsilons[j] >= data[i]){
          k =  j; break;
          }
        }
      vector[i] =Ys(k - 1) + (Ys(k) - Ys(k - 1)) * (data(i) - epsilons(k - 1)) / (epsilons(k) - epsilons(k - 1));
      }
    return vector;
    }

// [[Rcpp::export]]
NumericMatrix makeRatSplines(NumericMatrix matrice, NumericVector h) {
  int i, j, k, nr=matrice.nrow(), nc=matrice.ncol();
  NumericVector norm(nc);
  NumericMatrix matriceb(nr, nc);
  double s;
  for (j=k=0; j<nc; j++) {
    for (i=0, s=0.0; i<nr; i++,k++) {
      matriceb[k] = matrice[k]*h[i];
      s += matriceb[k];
      }
    norm[j] = s;
  }
  for (j=k=0; j<nc; j++) {
    for (i=0; i<nr; i++,k++) {
      matriceb[k] = matriceb[k]/norm[j];
    }
  }
  return matriceb;
}


// [[Rcpp::export]]
NumericVector makeWeights(NumericMatrix x){
  int i, nr = x.nrow();
  double r, s, t;
  NumericVector coseni(nr);
  for(i=0;i<nr-2;i++){
    r = (x(i,1)-x(i+1,1))*(x(i+2,1)-x(i+1,1))+(x(i,2)-x(i+1,2))*(x(i+2,2)-x(i+1,2));
    s = pow((x(i,1)-x(i+1,1)),2)+pow((x(i,2)-x(i+1,2)),2);
    t = pow((x(i+2,1)-x(i+1,1)),2)+pow((x(i+2,1)-x(i+1,2)),2);
    coseni(i) = r/(s*t);
    }
  return coseni;
}

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

// [[Rcpp::export]]
List findNewDimKnot(
    IntegerVector dcumFixedDim_Dim,
    NumericVector Dim_weights,
    NumericVector Dim_oldknots,
    NumericMatrix matrFixedDim,
    int Dim_index) 
{
  int u = dcumFixedDim_Dim.size(); // Number of clusters
  bool flagDim = false; // Flag for negative weights
  double Dim_newknot = NumericVector::get_na();  
  double weightDim = NumericVector::get_na();    
  
  int Dim_index_Cpp = Dim_index - 1;  // Convert R index (1-based) to C++ (0-based)
  int dcumInf = -1, dcumSup = -1;  // Initialize with invalid values
  
  if (Dim_index_Cpp < 0 || Dim_index_Cpp >= matrFixedDim.ncol()) {
    stop("Dim_index is out of bounds.");
  } 
  
  for (int i = 0; i < u; ++i) {
    if (is_true(all(Dim_weights < 0))) {
      flagDim = true;
      break;
    }
    
    // Find the index of the cluster with the highest weight
    int best_index = which_max(Dim_weights);
    
    if (best_index < 0 || best_index >= dcumFixedDim_Dim.size()) {
      stop("Invalid cluster index detected.");
    } 
    
    // Determine boundaries indexes
    dcumInf = (best_index == 0) ? 0 : dcumFixedDim_Dim[best_index - 1];
    dcumSup = dcumFixedDim_Dim[best_index] - 1;
    
    if (dcumInf >= matrFixedDim.nrow() || dcumSup >= matrFixedDim.nrow() || dcumInf > dcumSup) {
      stop("Indexing error: dcumInf or dcumSup out of range.");
    } 
    
    // Calculate superior and inferior Dim-bounds
    double sup = matrFixedDim(dcumSup, Dim_index_Cpp);
    double inf = matrFixedDim(dcumInf, Dim_index_Cpp);
    
    // Compute weighted average using std::inner_product and std::accumulate
    NumericVector weights = matrFixedDim(_, 2);
    NumericVector values = matrFixedDim(_, Dim_index_Cpp);
    
    Dim_newknot = std::inner_product(
      weights.begin() + dcumInf, weights.begin() + dcumSup + 1, 
      values.begin() + dcumInf, 
      0.0
    ) / std::accumulate(weights.begin() + dcumInf, weights.begin() + dcumSup + 1, 0.0);
    
    // Validate new knot conditions
    bool cond1 = (dcumSup - dcumInf) != 0;
    bool cond2 = !is_true(any((Dim_oldknots >= inf) & (Dim_oldknots <= sup)));
    bool cond3 = !is_true(any(abs(inf - Dim_oldknots) < 1e-12));
    
    if ((cond1 && cond2) || (!cond1 && cond3)) {
      weightDim = Dim_weights[best_index];
      break;
    } else {  
      Dim_weights[best_index] = R_NegInf;
    }
  }  
  
  return List::create(
    Named("Dim.newknot") = Dim_newknot,
    Named("weightDim") = weightDim,
    Named("flagDim") = flagDim,
    Named("dcumInf") = dcumInf + 1,  // Convert to 1-based indexing
    Named("dcumSup") = dcumSup + 1  // Convert to 1-based indexing
  );
} 




