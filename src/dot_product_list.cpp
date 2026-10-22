#include <Rcpp.h>
#include <math.h>
#include <cmath>
#include <norm_dot_product.cpp>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> normDotProductList(std::vector<double> peaks_A, std::vector<double> ints_A, 
                                       std::vector<std::vector<double> > peaks_B, 
                                       std::vector<std::vector<double> > ints_B, double bin_size) {
  int n_B = peaks_B.size();
  std::vector<double> similarities;
  
  // compare spectra A with all in B list
  for (int i = 0; i < n_B; i++) {
    similarities.push_back(normDotProduct(peaks_A, ints_A, peaks_B[i],ints_B[i], bin_size));
  }
  
  return (similarities);
}

// performs the normDotProductShif between spec A and the list of spec B
// only returns the similarity values >= sim_cutoff (default 0.1)
// return a list with 3 items:
// (0) the list of valid cosines (>=0.1), (1) the list with their number of matches 
// and (2) the list with their index+1 in spec B list (R format - idx starts in 1 and not in 0)
// The spectra are read directly from the R objects (no copy to std::vector) and 
// the scratch buffers of normDotProductShiftCore are reused for all the spec B, 
// so there is no memory allocation during the comparisons.
// peaks_B and ints_B are lists of numeric vectors, mzs_B is a numeric vector 
// [[Rcpp::export]]
std::vector<std::vector<double>> normDotProductShiftList(
    const NumericVector& peaks_A, const NumericVector& ints_A,
    double mz_A, const List& peaks_B, 
    const List& ints_B, const NumericVector& mzs_B,
    double bin_size, double max_shift,  double sim_cutoff = 0.1) {
  // double ISO_MASS = 1.0033;  // mass (13C) - mass (12C)
  const int n_B = peaks_B.size();
  if (ints_B.size() != n_B || mzs_B.size() != n_B) {
    stop("peaks_B, ints_B and mzs_B must have the same length (%d, %d and %d).", 
         n_B, (int) ints_B.size(), (int) mzs_B.size());
  }
  std::vector<double> similarities, matches, valid_sim_idx;
  ShiftScratch ws; // scratch buffers reused in all comparisons
  const double* peaks_A_ptr = peaks_A.begin();
  const double* ints_A_ptr = ints_A.begin();
  const int n_A = peaks_A.size();
 
  // compare spectra A with all in B list
  for (int i = 0; i < n_B; i++) {
    double mz_diff = mz_A-mzs_B[i];
    
    // ignore shifts greater than the max shift allowed, default to 200 Da
    if (std::fabs(mz_diff) > max_shift) {
      mz_diff = 0.0;
    }
    // retrieve similarity value - cosine - and number of matched peaks 
    // (the vectors of the list are not copied)
    const NumericVector peaks_Bi = peaks_B[i];
    const NumericVector ints_Bi = ints_B[i];
    double cosine;
    int matched_peaks;
    normDotProductShiftCore(peaks_A_ptr, ints_A_ptr, n_A,
                            peaks_Bi.begin(), ints_Bi.begin(), (int) peaks_Bi.size(),
                            bin_size, mz_diff, ws, cosine, matched_peaks);
    // check if the similarity value is greater or equal than the cut of 0.1
    // if yes, return it and add the current spec B idx to the valid list, 
    // otherwise ignore it - which is equal to setting it to zero
    if (cosine >= sim_cutoff) {
      similarities.push_back(cosine);
      matches.push_back(matched_peaks);
      valid_sim_idx.push_back(i+1);
    }
    // otherwise dont return the sim and the matches:
  }
  // return similarities and matched peaks
  std::vector<std::vector<double> > out_sim_matches;
  out_sim_matches.push_back(similarities);
  out_sim_matches.push_back(matches);
  out_sim_matches.push_back(valid_sim_idx);
  return (out_sim_matches);
}
