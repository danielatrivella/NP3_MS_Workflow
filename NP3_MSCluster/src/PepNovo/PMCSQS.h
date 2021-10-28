#ifndef __PMCSQS_H__
#define __PMCSQS_H__

/*! \file PMCSQS.h

\brief	Contains the models for spectrum quality score (SQS) and parent mass correction PMC

\todo split these
*/

#include "RankBoost.h"
#include "PeakList.h"
#include "../MLlib/mlmodel.h"
#include "../MLlib/MlRankBoost.h"

#define MIN_SPECTRA_FOR_PMCSQS_MODEL 100



// The 3 SQS model makes binary decisions: good charge c spectrum or not good for charges 
// c=1,2,3.
// The PMC model assumes that we have a good spectrum and in the correct charge.
// 


/*! \struct PMCRankStats

\brief Holds statistics about peak mathces (pairs, pairs with neutral losses, etc.) for a given
precursor m/z. This struct is used for creating PMC and SQS features.
*/
struct PMCRankStats {

	void clear();

	PMCRankStats() { clear(); }
	PMCRankStats(const PMCRankStats& other);

	mass_t m_over_z;

	score_t rank_score;

	float numFragmentPairs;
	float numStrongFragmentPairs;
	float numCharge2FragmentPairs;
	float numStrongCharge2FragmentPairs;
	float numH2oLossFragmentPairs;
	float numCharge2H2oLossFragmentPairs;

	float intensityInFragmentPairs;
	float intensityInStrongPairs;
	float intensityInCharge2Pairs;
	float intensityInStrongCharge2Pairs;
	float intensityInH2oLossPairs;
	float intensityInCharge2H2oLossPairs;

	float meanOffsetPairs;
	float meanOffsetStrongPairs;
	float meanOffsetCharge2Pairs;
	float meanOffsetCharge2StrongPairs;
	float meanOffsetH2oPairs;
	float meanOffsetCharge2H2oPairs;

	bool ind_pairs_with_min_tol;
	bool ind_strong_pairs_with_min_tol;
	bool ind_c2_pairs_with_min_tol;
	bool ind_c2_strong_pairs_with_min_tol;


	vector<float> offsetPairsOderedByIntensity;
	vector<float> strong_offsetPairsOderedByIntensity;
	vector<float> charge2offsetPairsOderedByIntensity;

	float numPairsIso0, intensityInPairsIso0;
	float numPairsIso1, intensityInPairsIso1;
	float numPairsIso2, intensityInPairsIso2;

	float numCharge2PairsIso0, intensityInCharge2PairsIso0;
	float numCharge2PairsIso1, intensityInCharge2PairsIso1;
	float numCharge2PairsIso2, intensityInCharge2PairsIso2;
};


// per charge
struct PmcSqsChargeRes {
	
	PmcSqsChargeRes() : mz1(-1.0), score1(NEG_INF), mz2(-1.0), score2(NEG_INF), min_comp_prob(-1.0), sqs_prob(0.0) {}

	mass_t mz1;
	float score1;
	mass_t mz2;
	float score2;

	float min_comp_prob; // the minimal probaiblity when comparing to other charges
	float sqs_prob;
};

// A struct that holds info for the dynamic programming table used in SQS
struct  DPColumn {
	short pointers[Val+1];
	short fPointers[Val+1];	
};

/*!
\class PMCSQS_Scorer

\brief holds all models for both SQS and PMC

\todo Split into two classes (currently they share some data structures so that needs to be worked out)
\todo Sdd operaotrs to features in models for crap vs. all good.
*/
class PMCSQS_Scorer {
public:

	PMCSQS_Scorer() : maximalChargeWithModels_(0),
				fragmentPairSumOffset(NEG_INF), bin_increment(NEG_INF), 
				indInitializedPmc_(false), indInitializedSqs_(false), config_(NULL), 
				currentSpectrumTotalIntensity_(0), currentSpectrumStrongIntensity_(0),
				currentSpectrumNumStrongPeaks_(0) {};

	
	bool initializeForCurrentSpectrum(const Config* config, const PeakList& pl);

	double calculateSqsScore(const Config* config, const PeakList& pl, size_t* charge);

	float computePmcsqsResultsForSpectrum(const Config* config,
										  const PeakList& pl,
										  vector<PmcSqsChargeRes>& results);

	void computeBestMzValuesForCharge(const PeakList& pl, 
										int charge,
										mass_t precursorMassTolerance,
										PmcSqsChargeRes& result);

	int getNumSizes() const { return sqsMassThresholds_.size() +1; }
	int getMaximalChargeWithModels() const { return maximalChargeWithModels_; }

	/*! \fn getSqsSizeIndex
	\brief computes the size index for sqs models according to the preursor m/z. 
	
	\note This size index is not the same one used in the other scoring models, 
	since it uses one set of thresholds for all charges (at the SQS stage we do not
	assume to know the precursor charge).

	@param precursorMz the precursor m/z written in the spectrum file
	*/
	int  getSqsSizeIndex(mass_t precursorMz) const
	{
		int i;
		for (i=0; i<sqsMassThresholds_.size(); i++)
			if (precursorMz<=sqsMassThresholds_[i])
				break;
		return i;
	}


	void   setFragmentPairSumOffset(mass_t offset) { fragmentPairSumOffset = offset; }
	mass_t get_frag_pair_sum_offset() const { return fragmentPairSumOffset; }

	void   setBinMassIncrement(mass_t inc) { bin_increment = inc; }
	mass_t get_bin_increment() const { return bin_increment; }


	const vector< vector<PMCRankStats> >& get_curr_spec_rank_pmc_tables() const { return currentSpectrumPmcTables_; }

	

	float computeBestMzAndChargeForSpectrum(const Config *config, const PeakList& pl,
							mass_t* mz1, int* charge1, float *prob1,
							mass_t* mz2, int* charge2, float *prob2,
							vector<PmcSqsChargeRes>* all_res = NULL);



	void selectPrecursorMassesAndCharges(const Config *config, 
										 const PeakList& pl,
										 vector<mass_t>& precursorMassesWith19,
										 vector<int>&    charges,
										 vector<PmcSqsChargeRes>* allResults=NULL);


	void test_pmc(Config *config, char *specs_file, int charge, mass_t min_mass=0,
		mass_t max_mass = POS_INF);

	

	void trainRegressionSqsModels(const Config* config, 
								   const SpectraAggregator& positiveSpectra,
								   const char* pathNegativeSpectraList,
								   double weightNegativeSpectra = 0.5);


	/*! \brief Trains the Precursor Mass Correction RankBoost models
	*/
	void trainPmcRankModels(const Config* config, 
							const SpectraAggregator& sa, 
							int specificCharge = 0, 
							bool overwrite=true);

	/*! \brief Reads the PMC models file (all models in the same file).
	*/
	bool readPmcModels(const Config *config, const char *file);



	bool read_pmc_rank_models(const Config *config, 
							  const char *file);

	void write_pmc_rank_models(const char *path) const;

	void write_sqs_models(const char *path) const;

	bool readSqsModels(const Config* config, const char* fileName);
	void writeSqsModels(const char* path) const;


	void benchmark_sqs(Config *config, char *list, char *anns);

	void benchmarkPm(const Config* config, const char* list);

	bool getIndInitializedPmc() const { return indInitializedPmc_; }
	bool getIndInitializedSqs() const { return indInitializedSqs_; }

	int get_max_model_charge() const { return maximalChargeWithModels_; }

/*	void output_filtered_spectra_to_mgfs(Config *config,
									 const vector<string>& files,
									 char *out_dir,
									 float filter_prob, 
									 int& num_written, 
									 int& num_read);*/

private:

	vector< vector<RankBoostModel * > >  pmc_rank_models; // charge / size thresholds

	vector< vector<MlModel*> >			 pmcRankBoostModels_; // charge / size thresholds
	vector< vector<mass_t> >             pmcMassThresholds_;
	vector< vector<mass_t> >			 pmcMzBiases_;


	vector<mass_t> sqsMassThresholds_; // thsese thresholds determine the model index of for SQS
									   // since we do not know the charge in advance, we cannot use
									   // the same size thresholds determined in the Config


	vector< MlModel* > regressionSqsModels_; // models for crap vis. good charge
											// model for charge 0 is crap against good from all charges
										
	vector<double> sqsRegressionNormalizingConstants_; // constants used to modify SQS scores so a score
													  // of 0.1 will keep 98% of the good spectra in each size


	vector< vector< vector<MlModel*> > > regressionChargeModels_; // models for charge1 vs. charge 2 (for slecting charge)
																  // dims [ size / charge1 / charge 2]


	int maximalChargeWithModels_; // the maximal charge for which we have a model

	mass_t fragmentPairSumOffset; // b+y or c+z - (PM+19)

	mass_t bin_increment;

	bool indInitializedPmc_;

	bool indInitializedSqs_;


	const Config* config_;

	float currentSpectrumTotalIntensity_;
	float currentSpectrumStrongIntensity_;
	int   currentSpectrumNumStrongPeaks_;

	vector<bool>  currentSpectrumStrongPeakFlags_;	   // indicator if the peak should be considered strong
	vector<bool>  currentSpectrumStrictIsotopeFlags_; // holds for each peak ind if there is any peak at -1
	vector<float> currentSpectrumIsotopeLevels_;

	vector< vector<PMCRankStats> > currentSpectrumPmcTables_;
	vector<PMCRankStats>		currentSpectrumBackgroundStats_;
	vector<PMCRankStats>		currentSpectrumMaximalValues_;
	
	void fillSqsDynamicProgrammingTable(const PeakList& pl, vector<DPColumn>& dp, int fargmentCharge ) const;

	void fillSqsSample(const PeakList& pl, MlSample& sample) const;

	/*! \brief fills the samples PMC samples for different masses.
	*/
	void fillRankboostPmcSamples(const PeakList& pl,
								 int charge,
								 vector<MlSample>& samples) const;

	void fillRankboostPmcSamplesOld(const PeakList& pl,
								 int charge,
								 vector<RankBoostSample>& samples) const;

	
	size_t findOptimalBinIdx(int true_mz_bin, int charge) const;


	/* \brief Selects which of the samples for the PMCR will be used in training.

	Takes a set of samples around the correct mass ([-3+5] every 0.1 Da.)
	Selects the bin of the correct mass as positive and a set from offseted m/z
	as negative samples. 
	*/
	void selectTrainingSampleIndexes(
		int charge,
		const vector<MlSample>& samples,
		const PeakList& pl,
		size_t& correctIndex,
		vector<size_t>& badPmcIndexes) const;

	size_t calcRankModelSizeIndex(int charge, mass_t pm_with_19) const
	{
		if (pmcMassThresholds_.size()<=charge)
		{
			cout << "Error: PMC does not support charge " << charge << endl;
			exit(1);
		}

		size_t i;
		for (i=0; i<pmcMassThresholds_[charge].size(); i++)
			if (pm_with_19<=pmcMassThresholds_[charge][i])
				break;
		return i;
	}


	// computes the sqs score using only the regression models 
	double calculateSqsScoreFromSample(const MlSample& sample, size_t sizeIndex) const;
	int    selectChargeForSample(const MlSample& sample, size_t sizeIndex,
								 double* p=0, int* runnerUp=0) const;

	void calculateCurrentSpectrumPmcValues(const PeakList& pl, mass_t binIncrement,
										   int selectedCharge = 0);

	void calculatePmcValuesForSqs(const PeakList& pl, mass_t binIncrement);

	void getSqsFeaturesFromPmcTabels(const PeakList& pl, 
									 vector< vector<float> >& sqs_featrues) const;


	void setPmcrMassThresholds(int option = 0); // option 1 , fixed for it data

	void set_sqs_mass_thresholds();


	// selects the class weights according to the input data
	// makes sure that the total weight assigned to the spectra from the negative spectra
	// (of all sizes) equals the weight given.
	/*! fn setClassWeightsAccordingToData
	\breif Sets the weights for the different paritions such that the total weight of
	all spectra is 1.0, and the total weight of the negative spectra is given as an input.

	@param positiveSpectra		 holds headers of all positive spectra
	@param negativeSpectra		 holds headers of all negative ("crap") spectra
	@param weightNegativeSpectra should be 0-1.0 (something like 0.3-0.6 is typical)
	@param classWeights			 the weight assignments to all spectra
	*/
	void setClassWeightsAccordingToData(const SpectraAggregator& positiveSpectra,
										const SpectraAggregator& negativeSpectra,
										double	weightNegativeSpectra,
										vector< vector<double> >& classWeights,
										bool verbose=true) const;

	
	/*! \fn findSqsNormalizingConstant
	\brief Computes the normalizing constant that brings the SQS score to a state
	in which 98% of the good spectra have a score >=0.1.
	
	@param sizeIndex SQS size index of the precursor m/z for the model.
	@param params  pointer to a the parameters involved in model learning (MlTrainingContainer)
	@param ratioPositiveMissed portion of good spectra allowed under 0.1 (default 0.02)
	@param verbose set true if you want it to make more noise

	@return the sqs score threshold for which \e ratioPositiveMissed of the positive spectra have a lower score
	*/
	double findSqsNormalizingConstant(size_t sizeIndex,
									  MlTrainingContainer* params, 
									  double ratioPositiveMissed = 0.02,
									  bool verbose = true);

	/*! \fn findSqsNormalizingConstantForOneAgainstAll
	\brief Computes the normalizing constant that brings the SQS score to a state
	in which 98% of the good spectra have a score >=0.1. This function is used for
	a model that pits good spectra from all charges against the bad spectra.
	*/
	void findSqsNormalizingConstantForOneAgainstAll(size_t sizeIndex,
													MlTrainingContainer* params, 
													double ratioPositiveMissed = 0.02,
													bool verbose = true);
};





void fillPmcRankStatistics(int charge,
						  const mass_t singleChargePairOffset, // the sum of b+y or c+z
						  mass_t minusRange, 
						  mass_t plusRange,
						  mass_t increment,
						  const Config* config,
						  const PeakList& pl,
						  const vector<bool>& strongPeakIndicators,
						  const vector<float>& isotopicLevels,
						  const vector<bool>& strictIsotopicIndicators,
						  vector<PMCRankStats>& pmcStatisticsVector);




void calculatePmcRankStatisticsForMass(const PeakList& pl, 
									   mass_t single_charge_pair_sum,
									   mass_t tolerance, 
									   const vector<float>& iso_levels,
									   const vector<bool>& strict_iso_inds,
									   const vector<bool>& strongPeakInds,
									   PMCRankStats& statistics);


void calculateBackgroundStatistics(const mass_t single_charge_pair_sum, // the sum of b+y or c+z
						  const Config *config,
						  const PeakList& pl,
						  const vector<bool>& strong_inds,
						  const vector<float>& iso_levels,
						  const vector<bool>& strict_iso_inds,
						  PMCRankStats& pmc_stats_total);





void create_training_files(Config *config);

void init_PMC_feature_names(vector<string>& names);




#endif

