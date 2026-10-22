#include <Rcpp.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double normDotProduct(std::vector<double> peaks_A, std::vector<double> ints_A, 
                      std::vector<double> peaks_B, std::vector<double> ints_B, double bin_size) {
  int n_A = peaks_A.size(), n_B = peaks_B.size();
  
  //Rcout << n_A;
  
  double sum_ints_A = 0.0, sum_ints_B = 0.0, top_sum = 0.0; // sum(sa^2); sum(sb^2); sum(sai*sbi);
  int idxA = 0, idxB = 0;
  
  while (idxA < n_A && idxB < n_B)
  {
    const double diff = peaks_A[idxA] - peaks_B[idxB];
    //Rcout << diff;
    
    if (diff <= 0) // pB > pA; increment pA
    {
      if (diff + bin_size >= 0)
      {
        sum_ints_A += ints_A[idxA] * ints_A[idxA];
        sum_ints_B += ints_B[idxB] * ints_B[idxB];
        top_sum += ints_A[idxA] * ints_B[idxB];
        idxA++;
        idxB++;
      } else {
        sum_ints_A += ints_A[idxA] * ints_A[idxA];
        idxA++;
      }
    } else {
      if (diff <= bin_size)
      {
        sum_ints_A += ints_A[idxA] * ints_A[idxA];
        sum_ints_B += ints_B[idxB] * ints_B[idxB];
        top_sum += ints_A[idxA] * ints_B[idxB];
        idxA++;
        idxB++;
      } else {
        sum_ints_B += ints_B[idxB] * ints_B[idxB];
        idxB++;
      }
    }
  }
  
  for (int i = idxA; i < n_A; i++) {
    sum_ints_A += ints_A[i] * ints_A[i];
  }
  for (int i = idxB; i < n_B; i++) {
    sum_ints_B += ints_B[i] * ints_B[i];
  }
  
  return (top_sum / sqrt(sum_ints_A * sum_ints_B));
}

// [[Rcpp::export]]
double normDotProductTrim(std::vector<double> peaks_A, std::vector<double> ints_A, 
                      std::vector<double> peaks_B, std::vector<double> ints_B, double bin_size, double trim_mz) {
  int n_A = peaks_A.size(), n_B = peaks_B.size();
  
  //Rcout << n_A;
  
  double sum_ints_A = 0.0, sum_ints_B = 0.0, top_sum = 0.0; // sum(sa^2); sum(sb^2); sum(sai*sbi);
  int idxA = 0, idxB = 0;
  
  while (idxA < n_A && idxB < n_B && peaks_A[idxA] < trim_mz && peaks_B[idxB] < trim_mz)
  {
    const double diff = peaks_A[idxA] - peaks_B[idxB];
    //Rcout << diff;
    
    if (diff <= 0) // pB > pA; increment pA
    {
      if (diff + bin_size >= 0)
      {
        sum_ints_A += ints_A[idxA] * ints_A[idxA];
        sum_ints_B += ints_B[idxB] * ints_B[idxB];
        top_sum += ints_A[idxA] * ints_B[idxB];
        idxA++;
        idxB++;
      } else {
        sum_ints_A += ints_A[idxA] * ints_A[idxA];
        idxA++;
      }
    } else {
      if (diff <= bin_size)
      {
        sum_ints_A += ints_A[idxA] * ints_A[idxA];
        sum_ints_B += ints_B[idxB] * ints_B[idxB];
        top_sum += ints_A[idxA] * ints_B[idxB];
        idxA++;
        idxB++;
      } else {
        sum_ints_B += ints_B[idxB] * ints_B[idxB];
        idxB++;
      }
    }
  }
  
  for (int i = idxA; i < n_A; i++) {
    if (peaks_A[i] > trim_mz)
      break;
    sum_ints_A += ints_A[i] * ints_A[i];
  }
  for (int i = idxB; i < n_B; i++) {
    if (peaks_B[i] > trim_mz)
      break;
    sum_ints_B += ints_B[i] * ints_B[i];
  }
  
  return (top_sum / sqrt(sum_ints_A * sum_ints_B));
}

// ---------------------------------------------------------------------------
// shifted cosine function
// Matches exactly equal peak masses first and for the remaining try to match with shifted peak masses
// mzShift: the exact mz difference between the precursor mass of the spectra A and B 
//   (not the absolute, because needs to guarantee that spectra A have a bigger precursor mass); 
//   recommended range is between 12 and 200, outside this range set it to 0
//
// The computation is done by normDotProductShiftCore, that works directly on the 
// data pointers (no vector copies) and reuses the scratch buffers between calls, 
// so it can be called many times (e.g. from normDotProductShiftList) without any 
// memory allocation. The arithmetic operations and their order are the same of the 
// previous implementation, the results are identical.
// Both spectra peaks must be sorted by increasing m/z.

// BEGIN_SHIFT_CORE
// scratch buffers reused between the calls
struct ShiftScratch {
  std::vector<int> idxA_nomatch;  // idxs of the peaks of A without an exact match
  std::vector<int> idxB_nomatch;  // idxs of the peaks of B without an exact match
  std::vector<char> B_matched;    // flags the idxB_nomatch entries matched with a shift
};

// out_cos: the shifted cosine value rounded to 3 decimal places
// out_matches: the number of common/matched peaks
inline void normDotProductShiftCore(const double* peaks_A, const double* ints_A, int n_A,
                                    const double* peaks_B, const double* ints_B, int n_B,
                                    double bin_size, double mzShift, ShiftScratch& ws,
                                    double& out_cos, int& out_matches) {
  // if precursor mass of A < precurosr mass of B -> mzShift is negative then shift spectra
  // to garantee that spectra A have a bigger precursor spectra
  // (swap the pointers instead of copying the vectors)
  if (mzShift < 0.0) {
    std::swap(peaks_A, peaks_B);
    std::swap(ints_A, ints_B);
    std::swap(n_A, n_B);
    mzShift = mzShift * -1;
  }
  
  std::vector<int>& idxA_nomatch = ws.idxA_nomatch;
  std::vector<int>& idxB_nomatch = ws.idxB_nomatch;
  idxA_nomatch.clear();
  idxB_nomatch.clear();
  
  int i, matched_peaks = 0;
  double sum_ints_A = 0.0, sum_ints_B = 0.0, top_sum = 0.0; // sum(sa^2); sum(sb^2); sum(sai*sbi);
  int idxA = 0, idxB = 0;
  
  // match equal peak masses
  while (idxA < n_A && idxB < n_B)
  {
    const double diff = peaks_A[idxA] - peaks_B[idxB];
    
    if (diff <= 0) // pB > pA; increment pA
    {
      if (diff + bin_size >= 0)
      {
        sum_ints_A += ints_A[idxA] * ints_A[idxA];
        sum_ints_B += ints_B[idxB] * ints_B[idxB];
        top_sum += ints_A[idxA] * ints_B[idxB];
        matched_peaks++;
        idxA++;
        idxB++;
      } else {
        idxA_nomatch.push_back(idxA);
        idxA++;
      }
    } else {
      if (diff <= bin_size)
      {
        sum_ints_A += ints_A[idxA] * ints_A[idxA];
        sum_ints_B += ints_B[idxB] * ints_B[idxB];
        top_sum += ints_A[idxA] * ints_B[idxB];
        matched_peaks++;
        idxA++;
        idxB++;
      } else {
        idxB_nomatch.push_back(idxB);
        idxB++;
      }
    }
  }
  
  // add no match indexes to the list
  for (i = idxA; i < n_A; i++) {
    idxA_nomatch.push_back(i);
  }
  const int nomatch_A = (int) idxA_nomatch.size();
  for (i = idxB; i < n_B; i++) {
    idxB_nomatch.push_back(i);
  }
  const int nomatch_B = (int) idxB_nomatch.size();
  
  // try to match the not matched masses with a mzShift, 
  // if it is greater than the bin size
  ws.B_matched.assign(nomatch_B, 0);
  if (mzShift > bin_size) {
    int n_B_left = nomatch_B; // number of not matched peaks of B 
    // first entry of idxB_nomatch that can still match: the entries before it are too small
    // (pA - pB > mzShift + bin_size) for the current peak of A and for all the next ones 
    // (the peaks of A are sorted), so they are never scanned again
    int j_start = 0;
    for (i = 0; i < nomatch_A && n_B_left > 0; i++) {
      idxA = idxA_nomatch[i];
      const double peak_A = peaks_A[idxA];
      int j = j_start;
      bool matched = false;
      
      // find a match on peaks_B for idxA
      if (mzShift > 0.0) {
        while (j < nomatch_B) {
          if (ws.B_matched[j]) { // already matched, not in the list anymore
            j++;
            j_start = j;
            continue;
          }
          idxB = idxB_nomatch[j];
          
          const double diff = peak_A - peaks_B[idxB];
          if (diff <= 0) // pB > pA; increment pA
          {
            // we are looking for shifted peaks in A, not in B. Then we need pA > pB
            break; // no match for idxA
          } else { // pA > pB;
            const double diff2 = std::fabs(diff - mzShift);  // pA -pB - mzShift
            if (diff2 <= bin_size) // match with mzShift
            {
              sum_ints_A += ints_A[idxA] * ints_A[idxA];
              sum_ints_B += ints_B[idxB] * ints_B[idxB];
              top_sum += ints_A[idxA] * ints_B[idxB];
              matched_peaks++;
              ws.B_matched[j] = 1; // remove idxB from the not matched list
              n_B_left--;
              matched = true;
              break;
            } else if (diff > mzShift) {
              j++;  // no match for idxB,  increment pB
              j_start = j;
            } else {  // diff < mzShift
              break; // no match for idxA,  increment pA
            }
          }
        }
      }
      if (!matched) {
        sum_ints_A += ints_A[idxA] * ints_A[idxA];
      }
    }
  } else {
    i = 0;
  }
  
  // add the no matched peaks
  for (;i < nomatch_A; i++) {
    sum_ints_A += ints_A[idxA_nomatch[i]] * ints_A[idxA_nomatch[i]];
  }
  for (i = 0; i < nomatch_B; i++) {
    if (!ws.B_matched[i]) {
      sum_ints_B += ints_B[idxB_nomatch[i]] * ints_B[idxB_nomatch[i]];
    }
  }
  
  // The multipliers to maintain the norm would cancel each other here 
  // eliminating the need to compute it
  if (top_sum > 0.0) {
    /// round similarity value to 3 decimal places
    out_cos = std::round((top_sum / std::sqrt(sum_ints_A * sum_ints_B)) * 1000.0) / 1000.0;
  } else {
    out_cos = top_sum;
  }
  out_matches = matched_peaks;
}
// END_SHIFT_CORE

// return a vector with two values: 
// (0) the shifted cosine value rounded to 3 decimal places and (1) the number of common/matched peaks
// [[Rcpp::export]]
std::vector<double> normDotProductShift(const std::vector<double>& peaks_A, const std::vector<double>& ints_A,
                                        const std::vector<double>& peaks_B, const std::vector<double>& ints_B,
                                        double bin_size, double mzShift) {
  ShiftScratch ws;
  double cosine;
  int matched_peaks;
  normDotProductShiftCore(peaks_A.data(), ints_A.data(), (int) peaks_A.size(),
                          peaks_B.data(), ints_B.data(), (int) peaks_B.size(),
                          bin_size, mzShift, ws, cosine, matched_peaks);
  // return cosine and matched peaks
  std::vector<double> out_cos_matches;
  out_cos_matches.push_back(cosine);
  out_cos_matches.push_back(matched_peaks);
  return out_cos_matches;
}
