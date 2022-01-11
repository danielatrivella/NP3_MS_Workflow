#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
//#include <pybind11/complex.h>
//#include <pybind11/functional.h>
//#include <pybind11/chrono.h>

//#include <algorithm>

#include <math.h>
using namespace std;

// compute the absolute difference
double diff_abs(double a, double b) {
  if (a >= b) {
    return a-b;
  } else {
    return b-a;
  }
}


// join isotopic peaks, keep only the smaller mz and sum the intensities
// only carbon isotopic pattern is considered and at most the presence of 2 C13 atoms
// bin_size is the tolerance in mz for peaks
std::vector<std::vector<double> > joinIsotopicPeaks(std::vector<double> peaks,
                                                    std::vector<double> ints,
                                                    double bin_size) {
  int n = peaks.size();
  double ISO_MASS = 1.0033;

  std::vector<int> peaks_rm;
  double totalIntensity = 0;

  for (int prev = 0; prev < n; prev++) {
    totalIntensity += ints[prev];
    // check if prev was already removed, if yes increment it to next position
    if (std::find(peaks_rm.begin(), peaks_rm.end(), prev) != peaks_rm.end()) {
      //Rcpp::Rcout << "skip prev\n";
      continue;
    }
    for (int i = prev+1; i < n; i++) {
      //Rcpp::Rcout << "prev " << prev << " i " << i << " n " << n <<" \n";
      // check if i was already removed, if yes skip to next position
      if (std::find(peaks_rm.begin(), peaks_rm.end(), i) != peaks_rm.end()) {
        //Rcpp::Rcout << "skip i\n";
        continue;
      }
      // if (diff_abs(peaks[i] - peaks[prev], ISO_MASS) <= bin_size ||
      //     diff_abs(peaks[i] - peaks[prev], 2*ISO_MASS) <= bin_size) {
      if ((diff_abs(peaks[i] - peaks[prev], ISO_MASS) <= bin_size ||
          diff_abs(peaks[i] - peaks[prev], 2*ISO_MASS) <= bin_size) &&
          (2/3*ints[prev] >= ints[i])) {
        //Rcpp::Rcout << "join \n";
        // join peaks intensities, keep the smaller mz
        ints[prev] = ints[prev] + ints[i];
        peaks_rm.push_back(i);
      } else if ((peaks[i]-peaks[prev]) > 2*ISO_MASS + bin_size) {
        //Rcpp::Rcout << "increment prev \n";
        break;
      }
    }
  }
  //Rcpp::Rcout << "end peaks to remove m = " << peaks_rm.size() << "\n";
  sort(peaks_rm.begin(), peaks_rm.end());
  std::vector<double> out_p, out_i;
  double totalIntensity_out = 0;

  int m = peaks_rm.size();
  //Rcpp::Rcout << "m " << m << " \n";
  for (int i = 0, j = 0; i < n; i++) {
    // skip the peaks that were joined from the peaks and ints list
    if (j >= m || peaks_rm[j] != i)
    {
      //Rcpp::Rcout << "joined \n";
      out_p.push_back(peaks[i]);
      out_i.push_back(ints[i]);
      totalIntensity_out += ints[i];
    } else {
      //Rcpp::Rcout << "j " << j << " peaks_rm[j] " << peaks_rm[j] << " \n";
      j++;
    }
  }

  // assert that the total intensity is still the same
  if (diff_abs(totalIntensity,totalIntensity_out) > 0.00001) {
    cout << "Original total intensity of the peaks: " << totalIntensity << "\n";
    cout << "Final total intensity after joining isotopic peaks: " << totalIntensity_out << "\n";
    //Rcpp::stop("Error joining isotopic peaks: the peaks' total intensity was not preserved!\n");
  }
  //Rcpp::Rcout << "out \n peaks size " << out_i.size() << "\n";
  std::vector<std::vector<double> > out;
  out.push_back(out_p);
  out.push_back(out_i);

  return out;
}


// bin_size : tolerance in m/z for two peaks be considered the same and be merged in to a single one; their intensity will be summed and their mz will be averaged using their intensities as weights
// trim_mz : the precursor mz to be trimmed - All peaks in a +/- 20 Da around precursor ion mass are deleted. This removes the residual precursor ion, which is frequently observed in MS/MS spectra acquired on qTOFs. if trim_mz = -1 disable the trimming
// scale_method : x = 0 - ln; x != 0 - power
// if join_isotopic_peaks equals 1 also run the joinIsotopicPeaks
// [[Rcpp::export]]
std::vector<std::vector<double> > joinAdjPeaksScalee(std::vector<double> peaks,
                                                     std::vector<double> ints,
                                                     double bin_size, double trim_mz,
                                                     double scale_method, int join_isotopic_peaks) {
  int n = peaks.size();
  std::vector<double> ints_joined;

  int prev = 0;
  double totalIntensity = ints[prev];

  for(int i = prev+1; i < n; ++i) {
    // if the precursor mass is informed, remove the peaks in a +/- 20 Da around precursor ion mass
    if (trim_mz != -1 && diff_abs(peaks[i], trim_mz) < 20) {
      continue;
    }
    totalIntensity += ints[i];
    if (peaks[i] - peaks[prev] <= bin_size) {
      // join peaks proportional to their intensity
      const double intensitySum = ints[prev] + ints[i];
      const double ratio = ints[prev]/intensitySum;
      const double newMz = peaks[prev]*ratio + peaks[i]*(1-ratio);

      peaks[prev] = newMz;
      ints[prev] = intensitySum;
    } else {
      ints_joined.push_back(ints[prev]);
      prev++;
      peaks[prev] = peaks[i];
      ints[prev] = ints[i];
    }
  }
  ints_joined.push_back(ints[prev++]);

  // obtain min intensity to filter very low intensity peaks (without a window)
  // before joining isotopic peaks
  // do not filter spectra with few peaks
  double min_intensity;
  if (ints_joined.size() > 10) {
    sort(ints_joined.begin(), ints_joined.end());
    min_intensity = ints_joined[2*prev/3]*0.001;
  } else {
    min_intensity = 0.0;
  }

  peaks.erase(peaks.begin()+prev,peaks.end());
  ints.erase(ints.begin()+prev,ints.end());

  // join isotopic peaks
  if (join_isotopic_peaks == 1) {
    std::vector<std::vector<double> > out_iso;
    out_iso = joinIsotopicPeaks(peaks, ints, bin_size);
    peaks = out_iso[0];
    ints = out_iso[1];
  }

  // store the joined peaks concatenated with the intensities
  std::vector<double> out_p, out_i;

  // rmv peaks with intensities bellow the min
  for (int i = 0; i < peaks.size(); i++)
  {
    if (ints[i] > min_intensity)
    {
      out_p.push_back(peaks[i]);
      out_i.push_back(ints[i]);
    }
    else  // ints[i] > min_intensity
    {
      totalIntensity = totalIntensity - ints[i];
    }
  }

  std::vector<std::vector<double> > out;
  // scale the intensity values
  double multVal = 1000.0 / totalIntensity;
  for (int i = 0; i < out_i.size(); i++) {
    // compute the scaled intensity
    if (scale_method == 0.0)
    {
      out_i[i] = 1.0 + log(1.0 + multVal * out_i[i]); // ln scale
    } else {
      out_i[i] = pow(multVal * out_i[i], scale_method); // pow  scale
    }
  }
  out.push_back(out_p);
  out.push_back(out_i);

  return out;
}

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


// Match exactly peak masses first and for the remaining try to match with shifted peak masses
// mzShift: the exact mz difference beteween the precursor mass of the spectra A and B (not the absolute, because needs to garantee that spectra A have a bigger precursor mass); recommended range is between 12 and 100, outside this range set it to 0
// [[Rcpp::export]]
double normDotProductShift(std::vector<double> peaks_A, std::vector<double> ints_A,
                      std::vector<double> peaks_B, std::vector<double> ints_B,
                      double bin_size, double mzShift) {
  // if precursor mass of A < precurosr mass of B -> mzShift is negative then shift spectra
  // to garantee that spectra A have a bigger precursor spectra
  if (mzShift < 0.0) {
    peaks_A.swap(peaks_B);
    ints_A.swap(ints_B);
    mzShift = mzShift * -1;
  }

  int n_A = peaks_A.size(), n_B = peaks_B.size(), i;

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

  // Rcout << "nA " << n_A << "\n";
  // Rcout << "nB " << n_B << "\n";

  // try to match the not matched masses with a mzShift,
  // if it is greater than the bin size
  if (mzShift > bin_size) {
    for (i = 0; i < n_A && n_B > 0; i++) {
      idxA = idxA_nomatch[i];
      int j = 0, matched = 0;

      // find a match on peaks_B for idxA
      while (j < n_B && mzShift > 0.0) {
        idxB = idxB_nomatch[j];

        const double diff = peaks_A[idxA]-peaks_B[idxB];
        if (diff <= 0) // pB > pA; increment pA
        {
          // we are looking for shifted peaks in A, not in B. Then we need pA > pB
          break; // no match for idxA
          // diff = abs(diff + mzShift);
          // if (diff <= bin_size) // match with mzShift
          // {
          //   sum_ints_A += ints_A[idxA] * ints_A[idxA];
          //   sum_ints_B += ints_B[idxB] * ints_B[idxB];
          //   top_sum += ints_A[idxA] * ints_B[idxB];
          //   idxB_nomatch.erase(idxB_nomatch.begin() + j);
          //   n_B--;
          //   matched = 1;
          //   break;
          // } else if (diff > mzShift) {
          //   break; // no match for idxA
          // } else {
          //   j++;  // no match for idxB
          // }
        } else { // pA > pB;
          const double diff2 = abs(diff - mzShift);  // pA -pB - mzShift
          if (diff2 <= bin_size) // match with mzShift
          {
            sum_ints_A += ints_A[idxA] * ints_A[idxA];
            sum_ints_B += ints_B[idxB] * ints_B[idxB];
            top_sum += ints_A[idxA] * ints_B[idxB];
            idxB_nomatch.erase(idxB_nomatch.begin() + j);
            n_B--;
            matched = 1;
            break;
          } else if (diff > mzShift) {
            j++;  // no match for idxB,  increment pB
          } else {  // diff < mzShift
            break; // no match for idxA,  increment pA
          }
          // } else  {
          //   j++;  // no match for idxB,  increment pB
          // }
        }
      }
      if (matched == 0) {
        sum_ints_A += ints_A[idxA] * ints_A[idxA];
      }
    }
  } else {
    i = 0;
  }

  // TODO heuristica para ordenar e uar as menos intensas se n_B for menor que N_A, se nao nao precisa pq vai usar tudo de todo jeito
  // add the no matched peaks
  for (;i < n_A; i++) {
    sum_ints_A += ints_A[idxA_nomatch[i]] * ints_A[idxA_nomatch[i]];
  }

  for (i = 0; i < n_B; i++) {
    sum_ints_B += ints_B[idxB_nomatch[i]] * ints_B[idxB_nomatch[i]];
  }

  // The multipliers to maintain the norm would cancel each other here eliminating the need to compute it
  if (top_sum > 0.0) {
    return (top_sum / sqrt(sum_ints_A * sum_ints_B));
  } else {
    return top_sum;
  }
}

PYBIND11_MODULE(dotprod, m) {
    m.doc() = "C++ written dotproduct calc function"; // optional module docstring

    m.def("normDotProduct", &normDotProduct, "Returs the dot product for given peaksA, peaksB and bin_size");

    m.def("normDotProductTrim", &normDotProductTrim, "Returs the dot product for peaksA, peaksB and bin_size with a mz_given trim");

    m.def("normDotProductShift", &normDotProductShift, "Returs the dot product for given peaksA, peaksB, bin_size and mzShift by precursor diff");

    m.def("joinAdjPeaksScale", &joinAdjPeaksScalee, "Join Adjacent peaks of a peak list");
}
