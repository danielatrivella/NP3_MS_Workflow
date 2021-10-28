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
std::vector<double> normDotProductShiftList(
    std::vector<double> peaks_A, std::vector<double> ints_A,
    double mz_A, std::vector<std::vector<double> > peaks_B, 
    std::vector<std::vector<double> > ints_B, std::vector<double> mzs_B,
    double bin_size) {
  // double ISO_MASS = 1.0033;  // mass (13C) - mass (12C)
  int n_B = peaks_B.size();
  std::vector<double> similarities;
  
  // compare spectra A with all in B list
  for (int i = 0; i < n_B; i++) {
    double mz_diff = mz_A-mzs_B[i];
    // double mz_diff = abs(mz_A-mzs_B[i]);
    // double mz_diff_iso = abs(abs(mz_diff) - ISO_MASS);
    // double mz_diff_iso2 = abs(abs(mz_diff) - 2*ISO_MASS);
    // //minimun diff in organic molecules is the carbon mass and max to small motifs <= 100Da
    // //except for isotopes, which the mz_diff equals ISO_MASS or 2*ISO_MASS
    // if (mz_diff_iso > bin_size && mz_diff_iso2 > bin_size &&
    //     (mz_diff < 12 - bin_size || mz_diff > 100 + bin_size)) {
    //   mz_diff = 0.0;
    // }
    similarities.push_back(normDotProductShift(peaks_A, ints_A,
                                               peaks_B[i],ints_B[i],
                                               bin_size, mz_diff));
  }
  
  return (similarities);
}