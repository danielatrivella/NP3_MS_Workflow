#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

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
    Rcpp::Rcout << "Original total intensity of the peaks: " << totalIntensity << "\n";
    Rcpp::Rcout << "Final total intensity after joining isotopic peaks: " << totalIntensity_out << "\n";
    Rcpp::stop("Error joining isotopic peaks: the peaks' total intensity was not preserved!\n");
  }
  //Rcpp::Rcout << "out \n peaks size " << out_i.size() << "\n";
  std::vector<std::vector<double> > out;
  out.push_back(out_p);
  out.push_back(out_i);
  
  return out;
}

// TODO add full sanity check after norm
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
    // if the precursor mass is informed, 
    // remove the peaks in a +/- 20 Da around precursor ion mass
    // without removing the precursor m/z
    if (trim_mz != -1 && diff_abs(peaks[i], trim_mz) < 20 && 
        diff_abs(peaks[i], trim_mz) > bin_size) {
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
  for (int i = 0; i < int(peaks.size()); i++)
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
  for (int i = 0; i < int(out_i.size()); i++) {
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

// bin_size peaks with mz diff < bin_size will be joined - to disable the join set bin_size to 0
// if trim_mz == -1 do not trim spectra; if trim_mz != -1 trim spectra by precursor mz +- 20
// scale_factor: 0 - ln; (0,] expoent of pow; and if 1 no scale; if -1 no norm
// join_isotopic_peaks == 1 join isotopic peaks of carbon C13 and C14, else do not join isotopic peaks 
// [[Rcpp::export]]
List readMgfPeaksList(std::string filePath, double bin_size, float trim_mz, 
                      double scale_factor, int join_isotopic_peaks) {
  FILE*  mgfStream = fopen(filePath.c_str(),"r");
  
  //double ISO_MASS = 1.0033;
  std::vector<double> mz, inty, mz_prec;
  std::vector<int> scans;
  std::vector< std::vector<double> > mzs, ints, spec_scale;
  
  int scanNumber = -1;
  double mOverZ = -1.0;
  
  char buffer[256];
  //assert(fileType_ == IFT_MGF);
  
  while (fgets(buffer, 256, mgfStream))
  {
    if (strncmp(buffer,"BEGIN IONS",10))
      continue;
    break;
  }
  
  long positionInFile = ftell(mgfStream);
  //float precursorIntensity_ = 0.0;
  //float firstPeakIntensity=0.0;
  
  if (positionInFile < 0) {
    Rcpp::stop("Bad skip position in mgf file! This can often be corrected by running unix2dos (or vice versa if appropriate)\n");
  }
  
  // read header info and first peak
  while (true)
  {
    if (! strncmp(buffer,"END IONS",8))
    {
      if (trim_mz != -1)
        trim_mz = mOverZ;
      
      //Rcpp::Rcout << "Scan number " << scanNumber << " \n";
      
      if (scale_factor != -1)
      {
        spec_scale = joinAdjPeaksScalee(mz, inty, bin_size, trim_mz, 
                                        scale_factor, join_isotopic_peaks);
        // append current peak list to the total list
        mzs.push_back(spec_scale[0]);
        ints.push_back(spec_scale[1]);
      } else if (join_isotopic_peaks == 1) {
        spec_scale = joinIsotopicPeaks(mz, inty, bin_size);
        // append current peak list to the total list
        mzs.push_back(spec_scale[0]);
        ints.push_back(spec_scale[1]);
      } else {
        // no join nor scale nor norm
        mzs.push_back(mz);
        ints.push_back(inty);
      }
      if (scanNumber == -1 || mOverZ == -1.0)
        Rcpp::stop("Error parsing the mgf file. The scan number or the pepmass from an ion is missing or was not properly read.");
      
      mz_prec.push_back(mOverZ);
      scans.push_back(scanNumber);
      // initialize vars
      mz.clear();
      inty.clear();
      scanNumber = -1;
      mOverZ = -1.0;
      
      // got to next ion
      while (fgets(buffer, 256, mgfStream))
      {
        if (strncmp(buffer,"BEGIN IONS",10))
          continue;
        break;
      }
    }
    // eof
    if( ! fgets(buffer, 256, mgfStream))
      break;
    
    //Rcpp::Rcout << buffer << " \n";
    
    if (!strncmp(buffer,"BEGIN IONS",10))
      continue;
    
    if (! strncmp(buffer,"SCAN=",5) )
    {
      //Rcpp::Rcout << "scan";
      if (sscanf(buffer+5,"%d",&scanNumber) != 1)
      {
        Rcpp::stop("Error: couldn't read scan number from mgf file!\n");
      }
      continue;
    }
    else if (! strncmp(buffer,"SCANS=",6) ) // this is the offical MGF field, only the first number is kept
    {
      //Rcpp::Rcout << "scans";
      if (sscanf(buffer+6,"%d",&scanNumber) != 1)
      {
        Rcpp::stop("Error: couldn't read scan number from mgf file!\n");
      }
      continue;
    }
    else if (! strncmp(buffer,"PEPMASS=",8))
    {
      if (sscanf(buffer+8,"%lf",&mOverZ) != 1)
      {
        Rcpp::stop("Error: couldn't read pepmass!");
      }
      continue;
    }
    else // is this a peak?
    {
      float mass = -1.0;
      float intensity = -1.0;
      
      // read all the peaks in the fragmentation list
      do
      {
        //Rcpp::Rcout << "peaks ? " << buffer << " \n";
        if (! strncmp(buffer,"END IONS",8) )
          break;
        
        std::istringstream is(buffer);
        is >> mass >> intensity;
        
        if (mass >0.0 && intensity>0.0)
        {
          mz.push_back(mass);
          inty.push_back(intensity);
        } else {  // this is probably a not expected fiel and not a peak, skip line
          break;
        }
      } while ( fgets(buffer, 256, mgfStream) );
    }
  }
  fclose(mgfStream);
  
  return List::create(_["SCANS"] = (scans), 
                      _["MZS"] =  wrap(mzs), 
                      _["INTS"] = wrap(ints),
                      _["PREC_MZ"] = (mz_prec));
}


// [[Rcpp::export]]
List readMgfHeader(std::string filePath) {
  FILE*  mgfStream = fopen(filePath.c_str(),"r");
  
  std::vector<double> mz_prec;
  std::vector<double> ms2_int;
  std::vector<double> base_peak_int;
  std::vector<double> peak_area, rt, rtmin, rtmax;
  std::vector<std::string> peak_id;
  std::vector<int> scans, num_peaks;
  
  int scanNumber = -1;
  double mOverZ = -1.0;
  std::string peakId_ = "";
  double retentionTime_ = -1.0, retentionTimeMax_ = -1.0, retentionTimeMin_ = -1.0;
  double peakArea_ = -1.0;
  double precursorIntensity_ = -1.0;
  double basePeakIntensity_ = -1.0;
  int numPeaks = -1;
  
  char buffer[256];
  //assert(fileType_ == IFT_MGF);
  
  while (fgets(buffer, 256, mgfStream))
  {
    if (strncmp(buffer,"BEGIN IONS",10))
      continue;
    break;
  }
  
  long positionInFile = ftell(mgfStream);
  //float precursorIntensity_ = 0.0;
  //float firstPeakIntensity=0.0;
  
  if (positionInFile < 0) {
    Rcpp::stop("Bad skip position in mgf file! This can often be corrected by running unix2dos (or vice versa if appropriate)\n");
  }
  
  // read header info and first peak
  while (true)
  {
    if (! strncmp(buffer,"END IONS",8))
    {
      // check mandatory fields
      if (scanNumber == -1 || mOverZ == -1.0)
        Rcpp::stop("Error parsing the mgf file. The scan number or the pepmass from an ion is missing or was not properly read.");
      
      scans.push_back(scanNumber);
      peak_id.push_back(peakId_);
      rt.push_back(retentionTime_);
      rtmin.push_back(retentionTimeMin_);
      rtmax.push_back(retentionTimeMax_);
      ms2_int.push_back(precursorIntensity_);
      mz_prec.push_back(mOverZ);
      peak_area.push_back(peakArea_);
      base_peak_int.push_back(basePeakIntensity_);
      num_peaks.push_back(numPeaks);
      
      // got to next ion
      while (fgets(buffer, 256, mgfStream))
      {
        if (strncmp(buffer,"BEGIN IONS",10))
          continue;
        break;
      }
      // initialize vars
      scanNumber = -1;
      mOverZ = -1.0;
      peakId_ = "";
      retentionTime_ = -1.0, retentionTimeMax_ = -1.0, retentionTimeMin_ = -1.0;
      peakArea_ = -1.0;
      precursorIntensity_ = -1.0;
      basePeakIntensity_ = -1.0;
      numPeaks = -1;
    }
    // eof
    if( ! fgets(buffer, 256, mgfStream))
      break;
    
    //Rcpp::Rcout << buffer << " \n";
    
    if (!strncmp(buffer,"BEGIN IONS",10))
      continue;
    
    if (! strncmp(buffer,"SCAN=",5) )
    {
      //Rcpp::Rcout << "scan";
      if (sscanf(buffer+5,"%d",&scanNumber) != 1)
      {
        Rcpp::stop("Error: couldn't read scan number from mgf file!\n");
      }
      continue;
    }
    else if (! strncmp(buffer,"SCANS=",6) ) // this is the offical MGF field, only the first number is kept
    {
      //Rcpp::Rcout << "scans";
      if (sscanf(buffer+6,"%d",&scanNumber) != 1)
      {
        Rcpp::stop("Error: couldn't read scan number from mgf file!\n");
      }
      continue;
    }
    else if (! strncmp(buffer,"NUM_PEAKS=",10) )
    {
      if (sscanf(buffer+10,"%d",&numPeaks) != 1)
      {
        Rcpp::stop("Error: couldn't read the number of peaks from mgf file!\n");
      }
      continue;
    }
    else if (! strncmp(buffer,"RT=",3) )
    {
      if (sscanf(buffer+3,"%lf",&retentionTime_) != 1)
      {
        Rcpp::stop("Error: couldn't read retention_time!\n");
      }
      continue;
    }
    else if (! strncmp(buffer,"RTINSECONDS=",12) ) // this is the official MGF field name
    {
      if (sscanf(buffer+12,"%lf",&retentionTime_) != 1)
      {
        Rcpp::stop("Error: couldn't read retention_time!");
      }
      continue;
    }
    // NP3 peak profile and area
    else if (! strncmp(buffer,"RTMIN=",6) ) // this is the NP3 MGF field name
    {
      if (sscanf(buffer+6,"%lf",&retentionTimeMin_) != 1)
      {
        Rcpp::stop("Error: couldn't read retention_time_min!");
      }
      continue;
    }
    else if (! strncmp(buffer,"RTMAX=",6) ) // this is the Np3 MGF field name
    {
      if (sscanf(buffer+6,"%lf",&retentionTimeMax_) != 1)
      {
        Rcpp::stop("Error: couldn't read retention_time_max!");
      }
      continue;
    }
    else if (! strncmp(buffer,"PEAK_AREA=",10) ) // this is the Np3 MGF field name
    {
      if (sscanf(buffer+10,"%lf",&peakArea_) != 1)
      {
        Rcpp::stop("Error: couldn't read peak_area!");
      }
      continue;
    }
    else if (! strncmp(buffer,"PEAK_ID=",8) ) // this is the Np3 MGF field name
    {
      int len = strlen(buffer)-1;
      if (buffer[len]=='\r' || buffer[len]=='\n' )
        buffer[len]='\0';
      if (buffer[len-1]=='\r' || buffer[len-1]=='\n' )
        buffer[len-1]='\0';
      peakId_ = buffer + 8;
      // keep reading the field while the end of line is not reached
      while (buffer[len]!='\0' && fgets(buffer, 256, mgfStream)) 
      {
        len = strlen(buffer)-1;
        if (buffer[len]=='\r' || buffer[len]=='\n' )
          buffer[len]='\0';
        if (buffer[len-1]=='\r' || buffer[len-1]=='\n' )
          buffer[len-1]='\0';
        peakId_ = peakId_ + buffer;
      }
      continue;
    }
    else if (! strncmp(buffer,"PRECURSOR_INTENSITY=",20) )
    {
      if (sscanf(buffer+20,"%lf",&precursorIntensity_) != 1)
      {
        Rcpp::stop("Error: couldn't read cluster size!");
      }
      continue;
    }
    else if (! strncmp(buffer,"PEPMASS=",8))
    {
      if (sscanf(buffer+8,"%lf",&mOverZ) != 1)
      {
        Rcpp::stop("Error: couldn't read pepmass!");
      }
      continue;
    }
    else // is this a peak?
    {
      float mass = -1.0;
      float intensity = -1.0;
      
      // read all the peaks in the fragmentation list
      do
      {
        //Rcpp::Rcout << "peaks ? " << buffer << " \n";
        if (! strncmp(buffer,"END IONS",8) )
          break;
        
        std::istringstream is(buffer);
        is >> mass >> intensity;
        
        if (mass >0.0 && intensity>0.0)
        {
          if (intensity > basePeakIntensity_) {
            basePeakIntensity_ = intensity;
          }
        } else {  // this is probably a not expected field and not a peak, skip line
          break;
        }
      } while ( fgets(buffer, 256, mgfStream) );
    }
  }
  fclose(mgfStream);
  
  return DataFrame::create(_["scans"] = (scans), 
                      _["mz"] =  (mz_prec), 
                      _["int"] = (ms2_int),
                      _["rt"] = (rt),
                      _["rtmin"] = (rtmin),
                      _["rtmax"] = (rtmax),
                      _["peak_area"] = (peak_area),
                      _["peak_id"] = (peak_id),
                      _["base_peak_int"] = (base_peak_int),
                      _["num_peaks"] = (num_peaks),
                      _["stringsAsFactors"] = false);
}
