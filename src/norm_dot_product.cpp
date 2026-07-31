#include <Rcpp.h>
#include <math.h>
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

// shifted cosine function
// Matches exactly equal peak masses first and for the remaining try to match with shifted peak masses
// mzShift: the exact mz difference between the precursor mass of the spectra A and B 
//   (not the absolute, because needs to guarantee that spectra A have a bigger precursor mass); 
//   recommended range is between 12 and 200, outside this range set it to 0
// return a vector with two values: 
// (0) the shifted cosine value rounded to 3 decimal places and (1) the number of common/matched peaks
// [[Rcpp::export]]
std::vector<double> normDotProductShift(std::vector<double> peaks_A, std::vector<double> ints_A,
                      std::vector<double> peaks_B, std::vector<double> ints_B,
                      double bin_size, double mzShift) {
  // if precursor mass of A < precurosr mass of B -> mzShift is negative then shift spectra
  // to garantee that spectra A have a bigger precursor spectra
  if (mzShift < 0.0) {
    peaks_A.swap(peaks_B);
    ints_A.swap(ints_B);
    mzShift = mzShift * -1;
  }
  
  int n_A = peaks_A.size(), n_B = peaks_B.size(), i, matched_peaks = 0;
  
  //Rcout << n_A;
  std::vector<int> idxA_nomatch, idxB_nomatch;
  
  double sum_ints_A = 0.0, sum_ints_B = 0.0, top_sum = 0.0; // sum(sa^2); sum(sb^2); sum(sai*sbi);
  int idxA = 0, idxB = 0;
  
  // match equal peak masses
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
        matched_peaks++;
        idxA++;
        idxB++;
      } else {
        idxA_nomatch.push_back(idxA);
        //sum_ints_A += ints_A[idxA] * ints_A[idxA];
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
        //sum_ints_B += ints_B[idxB] * ints_B[idxB];
        idxB_nomatch.push_back(idxB);
        idxB++;
      }
    }
  }
  
  // add no match indexes to the list
  for (i = idxA; i < n_A; i++) {
    idxA_nomatch.push_back(i);
  }
  n_A = idxA_nomatch.size();
  for (i = idxB; i < n_B; i++) {
    idxB_nomatch.push_back(i);
  }
  n_B = idxB_nomatch.size();
  //Rcout << "nA " << n_A << "\n";
  //Rcout << "nB " << n_B << "\n";
  
  // try to match the not matched masses with a mzShift, 
  // if it is greater than the bin size
  if (mzShift > bin_size) {
    for (i = 0; i < n_A && n_B > 0; i++) {
      idxA = idxA_nomatch[i];
      //Rcout << "idxA" << idxA << "\n";
      int j = 0, matched = 0;
      
      // find a match on peaks_B for idxA
      while (j < n_B && mzShift > 0.0) {
        idxB = idxB_nomatch[j];
        //Rcout << "idxB" << idxB << "\n";
        
        const double diff = peaks_A[idxA]-peaks_B[idxB];
        //Rcout << "peaks_A[idxA]" << peaks_A[idxA] << "\n";
        //Rcout << "peaks_B[idxB]" << peaks_B[idxB]<< "\n";
        if (diff <= 0) // pB > pA; increment pA
        {
          //Rcout << "no match" << "\n";
          // we are looking for shifted peaks in A, not in B. Then we need pA > pB
          break; // no match for idxA
        } else { // pA > pB;
          const double diff2 = abs(diff - mzShift);  // pA -pB - mzShift
          if (diff2 <= bin_size) // match with mzShift
          {
            //Rcout << "match!" << "\n";
            sum_ints_A += ints_A[idxA] * ints_A[idxA];
            sum_ints_B += ints_B[idxB] * ints_B[idxB];
            top_sum += ints_A[idxA] * ints_B[idxB];
            matched_peaks++;
            idxB_nomatch.erase(idxB_nomatch.begin() + j);
            n_B--;
            matched = 1;
            break;
          } else if (diff > mzShift) {
            //Rcout << "no match - increment pB" << "\n";
            j++;  // no match for idxB,  increment pB
          } else {  // diff < mzShift
            //Rcout << "no match - increment pA" << "\n";
            break; // no match for idxA,  increment pA
          }
        }
      }
      if (matched == 0) {
        //Rcout << "no match sum ints A" << "\n";
        sum_ints_A += ints_A[idxA] * ints_A[idxA];
      }
    }
  } else {
    i = 0;
  }
  
  // TODO heuristica para ordenar e uar as menos intensas se n_B for menor que N_A, se nao nao precisa pq vai usar tudo de todo jeito
  // add the no matched peaks
  //Rcout << "i" << i << "\n";
  for (;i < n_A; i++) {
    //Rcout << "idxA_nomatch" << idxA_nomatch[i] << "\n";
    sum_ints_A += ints_A[idxA_nomatch[i]] * ints_A[idxA_nomatch[i]];
  }
  for (i = 0; i < n_B; i++) {
    //Rcout << "idxB_nomatch" << idxB_nomatch[i] << "\n";
    sum_ints_B += ints_B[idxB_nomatch[i]] * ints_B[idxB_nomatch[i]];
  }
  // return cosine and matched peaks
  std::vector<double> out_cos_matches;
 
  // The multipliers to maintain the norm would cancel each other here 
  // eliminating the need to compute it
  if (top_sum > 0.0) {
    /// round similarity value to 3 decimal places
    out_cos_matches.push_back(std::round((top_sum / sqrt(sum_ints_A * sum_ints_B)) * 1000.0) / 1000.0);
  } else {
    out_cos_matches.push_back(top_sum);
  }
  out_cos_matches.push_back(matched_peaks);
  return out_cos_matches;
}

