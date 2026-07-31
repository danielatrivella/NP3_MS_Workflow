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

// performs the normDotProductShif between spec A and the list of spec B
// only returns the similarity values >= 0.1
// return a list with 3 items:
// (0) the list of valid cosines (>=0.1), (1) the list with their number of matches 
// and (2) the list with their index+1 in spec B list (R format - idx starts in 1 and not in 0)
// [[Rcpp::export]]
std::vector<std::vector<double>> normDotProductShiftList(
    std::vector<double> peaks_A, std::vector<double> ints_A,
    double mz_A, std::vector<std::vector<double> > peaks_B, 
    std::vector<std::vector<double> > ints_B, std::vector<double> mzs_B,
    double bin_size, double max_shift) {
  // double ISO_MASS = 1.0033;  // mass (13C) - mass (12C)
  int n_B = peaks_B.size();
  std::vector<double> similarities, matches, out_normDotProductShift,valid_sim_idx;
 
  // compare spectra A with all in B list
  for (int i = 0; i < n_B; i++) {
    double mz_diff = mz_A-mzs_B[i];
    
    // ignore shifts greater than the max shift allowed, default to 200 Da
    if (abs(mz_diff) > max_shift) {
      mz_diff = 0.0;
    }
    // retrieve similarity value - cosine - and number of matched peaks 
    out_normDotProductShift = normDotProductShift(peaks_A, ints_A,
                                                   peaks_B[i],ints_B[i],
                                                    bin_size, mz_diff);
    // check if the similarity value (first idx) is greater or equal than the cut of 0.1
    // if yes, return it and add the current spec B idx to the valid list, 
    // otherwise ignore it - which is equal to setting it to zero
    if (out_normDotProductShift[0] >= 0.1) {
      similarities.push_back(out_normDotProductShift[0]);
      matches.push_back(out_normDotProductShift[1]);
      valid_sim_idx.push_back(i+1);
    }
    // otherwise dont return the sim and the matches:
    //out_normDotProductShift[0] = 0.0;
    //out_normDotProductShift[1] = 0;
  }
  // return similarities and matched peaks
  std::vector<std::vector<double> > out_sim_matches;
  out_sim_matches.push_back(similarities);
  out_sim_matches.push_back(matches);
  out_sim_matches.push_back(valid_sim_idx);
  return (out_sim_matches);
}