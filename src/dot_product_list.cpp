#include <Rcpp.h>
#include <math.h>
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

// [[Rcpp::export]]
std::vector<std::vector<double>> normDotProductShiftList(
    std::vector<double> peaks_A, std::vector<double> ints_A,
    double mz_A, std::vector<std::vector<double> > peaks_B, 
    std::vector<std::vector<double> > ints_B, std::vector<double> mzs_B,
    double bin_size, double max_shift) {
  // double ISO_MASS = 1.0033;  // mass (13C) - mass (12C)
  int n_B = peaks_B.size();
  std::vector<double> similarities, matches, out_normDotProductShift;
 
  // compare spectra A with all in B list
  for (int i = 0; i < n_B; i++) {
    double mz_diff = mz_A-mzs_B[i];
    
    // ignore shifts greatex than the max shift allowed, default to 200 Da
    if (abs(mz_diff) > max_shift) {
      mz_diff = 0.0;
    }
    // retrieve similarities and matched peaks in another variable
    out_normDotProductShift = normDotProductShift(peaks_A, ints_A,
                                                   peaks_B[i],ints_B[i],
                                                    bin_size, mz_diff);
    if (out_normDotProductShift[0] < 0.1) {
      out_normDotProductShift[0] = 0.0;
      out_normDotProductShift[1] = 0;
    }
      
    similarities.push_back(out_normDotProductShift[0]);
    matches.push_back(out_normDotProductShift[1]);
  }
  // return similarities and matched peaks
  std::vector<std::vector<double> > out_sim_matches;
  out_sim_matches.push_back(similarities);
  out_sim_matches.push_back(matches);
  return (out_sim_matches);
}