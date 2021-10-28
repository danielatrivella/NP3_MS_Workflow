#ifndef __SPECTRUM_H__
#define __SPECTRUM_H__


#include "Config.h"
#include "PeakList.h"
#include "SpectraAggregator.h"


/////////////////////////////////////////////////////////////////
/// This class holds the common characteristics of a spectrum
/// A spectrum is an entity whose charge is known!
class Spectrum : public PeakList {
public:
	Spectrum() : originalPmWith19_(-1.0),	mOverZ_(-1.),
				 correctedPmWith19_(-1.0),	correctedPmScore_(NEG_INF), 
			     secondaryPmWith19_(-1.0),	secondaryPmScore_(NEG_INF),
				 maximalPeakMassToConsider_(-1), minimalPeakMass_(-1), 
				 maximalPeakMass_(-1),		     charge_(-1),
				 scanNumber_(MIN_INT),		   retentionTime_(-1),
				 clusterSize_(1),		   sizeIndex_(-1),             
				 title_("") {};


	Spectrum(const PeakList& pl) { copyFromPeakList(pl); }


	void copyFromPeakList(const PeakList& pl);

	bool	readSpectrum(const SpectraAggregator& sa,	
						 const SingleSpectrumHeader* header, 
						 bool filterSpectrum=true);

	void initializeSpectrum(void *ssf, bool indFilterSpectrum=true);
	void initializeSpectrum();


	const vector<float>& getLogLocalRanks() const  { return logLocalRanks_; }
	const vector<float>& getLogIntensities() const { return logIntensities_; }

	PeakRange findPeaksInRange(mass_t minRange, mass_t maxRange) const;	
	int		findPeakWithMaxIntensity(mass_t mass, mass_t tolerance) const;
	void	findIsotopicEnvelope(int p_idx, vector<float>& iso_intens, mass_t iso_tolerance, int charge) const;


	const SingleSpectrumHeader* getHeader() const { return header_; }

	mass_t      get_true_mass() const         { return  peptide_.get_mass(); }
	mass_t      get_true_mass_with_19() const { return  peptide_.get_mass() + MASS_OHHH; }
	mass_t      get_org_pm_with_19() const { return  originalPmWith19_; }
	void        set_org_pm_with_19( mass_t pm) { originalPmWith19_ = pm; }

	mass_t	    get_m_over_z() const { return mOverZ_; }
	void        set_m_over_z(mass_t mz) { mOverZ_ = mz; }

	void        calc_original_pm_with_19_for_different_charges(vector<mass_t>& pms_with_19) const;
	mass_t      get_corrected_pm_with_19() const { return correctedPmWith19_; }
	mass_t      get_secondary_pm_with_19() const { return secondaryPmWith19_; }
	score_t     get_corrected_pm_score() const { return correctedPmScore_; }
	score_t     get_secondary_pm_score() const { return secondaryPmScore_; }

	void        set_corrected_pm_with_19(mass_t pm) { correctedPmWith19_ = pm; }
	void        set_secondary_pm_with_19(mass_t pm) { secondaryPmWith19_ = pm; }
	void        set_corrected_pm_score(score_t score) { correctedPmScore_ = score; }
	void        set_secondary_pm_score(score_t score) { secondaryPmScore_ = score; }

	mass_t		getMaximalPeakMassToConsider() const { return maximalPeakMassToConsider_; }
	int         get_size_idx() const { return sizeIndex_; }
	void        set_size_idx(int s) { sizeIndex_ = s; }
	int         getCharge() const { return charge_; }
	void        setCharge(int c) { charge_ = c; }
	int         getScanNumber() const { return scanNumber_; }
	int         getClusterSize() const { return clusterSize_; }
	const Config*  getConfig() const { return config_; }
	const string&  getTitle() const { return title_; }
	void setTitle(const string& fname) { title_ = fname; }
	const Peptide& getPeptide()   const { return peptide_; }
	void  set_peptide(const Peptide& p) { peptide_ = p; }
	void  set_peptide(const string& pep_string) { peptide_.parseFromString(config_,pep_string); }
	void  convert_peptide_ILQK() { peptide_.convert_ILQK(config_); }

	mass_t      get_min_peak_mass() const { return minimalPeakMass_; }
	mass_t      get_max_peak_mass() const { return maximalPeakMass_; }
	int         get_peak_rank(int peak_idx) const { return ranks_[peak_idx]; }
	float       get_peak_iso_level(int peak_idx) const { return isotopicLevels_[peak_idx]; }
	float       get_peak_log_intensity(int peak_idx) const { return logIntensities_[peak_idx]; }
	float		getPeakLogLocalRank(int peakIdx) const { return logLocalRanks_[peakIdx]; }
	float		getPeakLogRandomProbability(int peakIdx) const { return logRandomProbabilities_[peakIdx]; }


	const vector<int>& get_strong_peak_idxs() const { return strongPeakIndexes_; }

	
	void output_as_MGF(ostream& os) const;

	
	void init_spectrum(bool perform_filtering = true);

	void create_random_spec_with_1_good_peak();

	void print_spectrum(ostream& os = cout) const;

	void print_expected_by(ostream& os = cout) const;

	void print_expected_fragment_peaks(vector<string>& frag_labels, ostream& os = cout) const;

	bool check_m_over_z_and_sequence(ostream& os = cout);

	void initWithHeaderInfo();


protected:

	mass_t  originalPmWith19_;
	mass_t	mOverZ_;
	mass_t	correctedPmWith19_;
	mass_t	correctedPmScore_;
	mass_t	secondaryPmWith19_;
	mass_t	secondaryPmScore_;
	mass_t  maximalPeakMassToConsider_;   // the maximal mass that can be aniticipated for a peak
	mass_t  minimalPeakMass_, maximalPeakMass_;
	int charge_;
	int scanNumber_;        // if originated from mzXML
	float  retentionTime_;  // if originated from mzXML
	int clusterSize_;      // if originated from a cluster
	int sizeIndex_;
	string title_;
	Peptide peptide_;

	vector<int> indexArray_;  // used to quickly locate peaks

	vector<int> ranks_;		 /// amongst all peaks (including isotopic peaks). Highest rank is 1

	vector<float> isotopicLevels_;  /// if IsoLevel>0 for a peak, it means this might be isotopic peak.
					               /// The higher the value, the more likely this is an isotopic peak.
	vector<float> logIntensities_;

	vector<float> logLocalRanks_;

	vector<float> logRandomProbabilities_;

	vector<int> strongPeakIndexes_; // holds the indexes of peaks that should be considered strong
							     // these generally are the top peaks in the spectrum that are
							    // not isotopic peaks, and some peaks that are locally strong



	// functions

	
	void copyHeaderInformation();

	void init_index_array();     // creates the fast lookup array for the peaks
	void join_adjacent_peaks();  // joins peaks that are close to each other
	void filter_peaks();         // removes the locally weak peaks
	void normalize_intensities();
	void calc_ranks();
	void calc_log_local_ranks(); // calcs the logarithm of the local rank for each
	void calc_isotope_levels();  // calcs for each peak the level it resembles an isotopic peak
	void select_strong_peaks();
	void set_log_random_probs(); // calcs for each peak the probability of observing it at random (based on
								 // the neighbor's distribution

	Spectrum& operator= (const Spectrum& other);
};


ostream& operator << (ostream& os, const Spectrum& spec);

/*****************************************************************
Returns the indices of all peaks that are within the mass range
******************************************************************/
inline PeakRange  Spectrum::findPeaksInRange(mass_t minRange, mass_t maxRange) const
{
	PeakRange pr;
	const int maxPeakIdx = getNumPeaks()-1;

	if (maxRange > maximalPeakMassToConsider_+0.1)
		maxRange = maximalPeakMassToConsider_+0.1;

	if (minRange<0)
		minRange=0;

	if (maxRange <= minRange)
		return pr;

	int indexMin=indexArray_[(int)(minRange)];
	int indexMax=indexArray_[(int)(maxRange)];
		
	if (indexMax<maxPeakIdx)
		indexMax++;

	while (peaks_[indexMin].mass < minRange && indexMin<indexMax)
		indexMin++;

	while (indexMax>0 && peaks_[indexMax].mass > maxRange)
		indexMax--;

	if (peaks_[indexMin].mass > maxRange || peaks_[indexMax].mass<minRange)
		return pr;

	pr.num_peaks = indexMax - indexMin+1;
	pr.low_idx  = indexMin;
	pr.high_idx = indexMax;

	return pr;
}




/***********************************************************
returns the idx of the peak that has the highest intensity
in the range (-1 is returned if no peak is found)
Balances between intensity and proximity to the expected position:
A peak at the edge needs to be at least 2 times stronger than
a peak exactly at the middle to be chosen.
************************************************************/
inline int Spectrum::findPeakWithMaxIntensity(mass_t expectedMass, mass_t tolerance) const
{
	PeakRange pr = findPeaksInRange(expectedMass-tolerance, expectedMass+tolerance);
	if (pr.num_peaks==0)
		return -1;

	if (pr.num_peaks ==1)
		return pr.low_idx;

	const float doubleTolerance = tolerance*2;

	mass_t maxWeightedValue = (doubleTolerance - fabs(peaks_[pr.low_idx].mass - expectedMass)) *
							   peaks_[pr.low_idx].intensity;

	int maxPeakIdx=pr.low_idx;
	
	// find closest peak to exp_pos
	for (int i=1; i<pr.num_peaks; i++)
	{
		const int peakIdx = pr.low_idx + i;

		mass_t weightedValue = (doubleTolerance - fabs(peaks_[peakIdx].mass - 
								expectedMass)) * peaks_[peakIdx].intensity;

		if (weightedValue>maxWeightedValue)
		{
			maxWeightedValue=weightedValue;
			maxPeakIdx = peakIdx;
		}
		else
			break;
	}
	return maxPeakIdx;
}



/// returns intensities at positions -2,-1,+1,+2
inline void Spectrum::findIsotopicEnvelope(int peakIdx, 
										   vector<float>& isotopicIntensities, 
										   mass_t isotopicTolerance, 
										   int charge) const
{
	const mass_t oneProtonMz = MASS_PROTON/(float)charge;
	const mass_t twoProtonMz = oneProtonMz + oneProtonMz;
	const mass_t peakMassMinusTolerance =peaks_[peakIdx].mass - isotopicTolerance;
	const mass_t minMinus2Mass =  peakMassMinusTolerance - twoProtonMz;
	const mass_t minMinus1Mass =  peakMassMinusTolerance - oneProtonMz;
	const mass_t minPlus1Mass  =  peakMassMinusTolerance + oneProtonMz;
	const mass_t minPlus2Mass  =  peakMassMinusTolerance + twoProtonMz;


	isotopicIntensities.clear();
	isotopicIntensities.resize(4,NEG_INF);

	int p=peakIdx;
	while (p>=0 && peaks_[peakIdx].mass > minMinus2Mass)
		p--;

	p++;

	const mass_t& peakMass = peaks_[peakIdx].mass;

	if (fabs(peaks_[p].mass- peakMass + twoProtonMz) <= isotopicTolerance)
		isotopicIntensities[0] = logIntensities_[p];

	while (p<peakIdx && peaks_[p].mass<minMinus1Mass)
		p++;

	
	if (fabs(peaks_[p].mass - peakMass + oneProtonMz) <= isotopicTolerance)
		isotopicIntensities[1] = logIntensities_[p];

	p=peakIdx+1;

	while (p < numPeaks_ && peaks_[p].mass<minPlus1Mass)
		p++;

	if (p < numPeaks_ && fabs(peaks_[p].mass - peakMass - oneProtonMz) <= isotopicTolerance)
		isotopicIntensities[2] = logIntensities_[p];

	while (p < numPeaks_ && peaks_[p].mass < minPlus2Mass)
		p++;

	if (p < numPeaks_ && fabs(peaks_[p].mass - peakMass - twoProtonMz) <= isotopicTolerance)
		isotopicIntensities[3] = logIntensities_[p];
}


#endif




