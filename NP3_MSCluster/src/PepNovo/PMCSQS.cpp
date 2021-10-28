#include "PMCSQS.h"
#include "PepNovo_auxfun.h"
#include "ME_REG.h"


/********************************************************************
Select a set of mzs nd charges that have high probabilities.
If unless the use_spectrum_mz is set to 1 in the config, it adds -1,+1  
(if there are too many mzs, not all these offfsets will be added)
Uses emprically set probability threhsolds to chose corrected pms.
*********************************************************************/
void PMCSQS_Scorer::selectPrecursorMassesAndCharges(
								const Config *config, 
								const PeakList& pl,
								vector<mass_t>& precursorMassesWith19,
								vector<int>&    charges,
								vector<PmcSqsChargeRes>* allResults)
{
	static const mass_t offsets[]={-MASS_ISO,MASS_ISO};
	static const int numAllOffsets = sizeof(offsets)/sizeof(mass_t);
	static const float minimalCompareProbForAdding  = 0.04;
	static const float minimalSqsProbForAdding		= 0.05;
	const mass_t halfTolerance = config->getTolerance() * 0.5;
	
	const SingleSpectrumHeader* const header = pl.getHeader();
	int		spectrumCharge	= header->getCharge();
	mass_t	spectrumMOverZ	= header->getMOverZ();
	const int originalSpectrumCharge = spectrumCharge;
	const mass_t originalSpectrumMOverZ = spectrumMOverZ;
	
	precursorMassesWith19.clear();
	charges.clear();

	int    specificCharge=0;
	mass_t specificMOverZ=0;
	
	if (config->get_use_spectrum_charge())
		specificCharge=spectrumCharge;

	static vector<PmcSqsChargeRes> results;
	int maximalprobabilityCharge=0;
	float maximalProbability = 0;

	if (spectrumCharge>0 && 
		config->get_use_spectrum_charge() && 
		config->get_use_spectrum_mz())
	{
		precursorMassesWith19.push_back( spectrumMOverZ * spectrumCharge - 
										(spectrumCharge-1)*MASS_PROTON );
		charges.push_back(spectrumCharge);
		maximalprobabilityCharge=spectrumCharge;
	}
	else
	{
		if (! this->indInitializedPmc_)
		{
			cout << "Error: no PMC model was read!" << endl;
			if (header->getCharge() <=0)
				cout << "Spectrum was read with charge <=0, so charge selection is needed!" << endl;
			exit(1);
		}

		if (config->get_use_spectrum_charge() && spectrumCharge >= this->pmc_rank_models.size())
		{
			// adjust charge and m/z to lowest charge modeled
			spectrumCharge = pmc_rank_models.size()-1;
			mass_t mass_with_19 = maximalprobabilityCharge * spectrumMOverZ - (maximalprobabilityCharge-1)*MASS_PROTON;
			spectrumMOverZ = (mass_with_19 + (spectrumCharge -1 ) * MASS_PROTON) / spectrumCharge;

			// TODO there should be a better way to do this correction
			SingleSpectrumHeader *nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
			nonConstHeader->setCharge(spectrumCharge);
			nonConstHeader->setMOverZ(spectrumMOverZ);
		}
		

	//	get_pmcsqs_results_for_spectrum(config,bs,res);
		computePmcsqsResultsForSpectrum(config, pl, results);

		if (specificCharge>0)
		{
			maximalprobabilityCharge = spectrumCharge;
			maximalProbability = results[spectrumCharge].min_comp_prob;
		}
		else
		{
			int c;
			for (c=1; c<results.size(); c++)
			{
				if (results[c].min_comp_prob>maximalProbability)
				{
					maximalprobabilityCharge = c;
					maximalProbability = results[c].min_comp_prob;
				}
			}
		}

		precursorMassesWith19.push_back(results[maximalprobabilityCharge].mz1 * 
										maximalprobabilityCharge - (maximalprobabilityCharge-1)*MASS_PROTON);
		charges.push_back(maximalprobabilityCharge);
		if (results[maximalprobabilityCharge].mz2>0)
		{
			precursorMassesWith19.push_back(results[maximalprobabilityCharge].mz2 * 
										maximalprobabilityCharge - (maximalprobabilityCharge-1)*MASS_PROTON);
			charges.push_back(maximalprobabilityCharge);
		}

	}

	// only added to the first pm_with_19, check that it doesn't overlap with others
	if (! config->get_use_spectrum_mz())
	{
		int numOffsets = numAllOffsets;
		int i;
		for (i=0; i<numOffsets; i++)
		{
			const mass_t pm_with_19 = precursorMassesWith19[0] + offsets[i];
			int j;
			for (j=0; j<precursorMassesWith19.size(); j++)
				if (charges[j]==maximalprobabilityCharge && 
					fabs(pm_with_19 - precursorMassesWith19[j]) < halfTolerance)
					break;
	
			if (j==precursorMassesWith19.size())
			{
				precursorMassesWith19.push_back(pm_with_19);
				charges.push_back(maximalprobabilityCharge);
			}
		}

		// add tol to -1 of 2nd
		if (charges[0]>=2)
		{
			const mass_t pm_with_19 = precursorMassesWith19[1] - MASS_PROTON;
			int j;
			for (j=0; j<precursorMassesWith19.size(); j++)
				if (charges[j]==maximalprobabilityCharge &&
					fabs(pm_with_19 - precursorMassesWith19[j]) < halfTolerance)
						break;
			
			if (j==precursorMassesWith19.size())
			{
				precursorMassesWith19.push_back(pm_with_19);
				charges.push_back(maximalprobabilityCharge);
			}
		}
		

		if (charges[0]>=3)
		{
			const mass_t pm_with_19 = precursorMassesWith19[1] + MASS_PROTON;
			int j;
			for (j=0; j<precursorMassesWith19.size(); j++)
				if (charges[j]==maximalprobabilityCharge &&
					fabs(pm_with_19 - precursorMassesWith19[j]) < halfTolerance)
						break;

			if (j==precursorMassesWith19.size())
			{
				precursorMassesWith19.push_back(pm_with_19);
				charges.push_back(maximalprobabilityCharge);
			}
		}
	}

	// add other charges if their comp probability and sqs probs are high enough
	if (specificCharge==0)
	{
		// find best charge
		int c; 
		float max_prob=-1.0;

		for (c=1; c<results.size(); c++)
		{
			if (c==maximalprobabilityCharge)
				continue;

			if (results[c].min_comp_prob > minimalCompareProbForAdding &&
				results[c].sqs_prob > minimalSqsProbForAdding)
			{
				precursorMassesWith19.push_back(results[c].mz1 * c - (c-1)*MASS_PROTON);
				charges.push_back(c);
			}
		}
	}

	if (allResults)
		*allResults = results;

	if (spectrumCharge != maximalprobabilityCharge)
	{
		// TODO there should be a better way to do this correction
		SingleSpectrumHeader *nonConstHeader = const_cast<SingleSpectrumHeader*>(header);
		nonConstHeader->setCharge(maximalprobabilityCharge);
		nonConstHeader->setMOverZ(originalSpectrumMOverZ);
	}
}


/***************************************************************************
Gives detailed m/z and prob values for the different charges.
Returns the highest sqs probability found for a spectrum.
****************************************************************************/
float PMCSQS_Scorer::computePmcsqsResultsForSpectrum(const Config* config,
													 const PeakList& pl,
													 vector<PmcSqsChargeRes>& results)
{
	ME_Regression_Sample sqs_sam;

	if (! initializeForCurrentSpectrum(config, pl))
		return 0; // corrupt file

//	calculate_curr_spec_pmc_values(bs,bin_increment);
	calculateCurrentSpectrumPmcValues(pl, bin_increment);

	const int sqs_size_idx = getSqsSizeIndex(pl.getHeader()->getMOverZ());
	if (results.size() <= maximalChargeWithModels_)
		results.resize(maximalChargeWithModels_+1);

	vector<float> minComparisonProbs(maximalChargeWithModels_+1,2.0);	// the minimal probability when used in a comparison model
	minComparisonProbs[0]=0.0;

/*	if (indInitializedSqs_)
	{
		const int maxSqsCharge = ( sqs_models.size()<maximalChargeWithModels_ ?
			sqs_models.size() : maximalChargeWithModels_);

		MlSample sample;
		fillSqsSample(pl, sample);
	
		int charge;
		for (charge=1; charge<=maxSqsCharge; charge++)
		{
			if (sqs_models[charge].size()<1 ||
				sqs_models[charge][0].size() <= sqs_size_idx)
				continue;
			
			float prob = sqs_models[charge][0][sqs_size_idx]

		//	cout << "SQS: " << charge << "\t" << prob;
			// correct prob
			prob += sqs_correction_factors[charge][sqs_size_idx];
			prob *= sqs_mult_factors[charge][sqs_size_idx];
			if (prob<0.0)
				prob=0.0;

			results[charge].sqs_prob=prob;

		//	 cout << "\t" << prob << endl;
		}
		
		int c;
		for (c=1; c<maxSqsCharge; c++)
		{
			if (sqs_models[c].size()<=0)
				continue;

			float comp_prob = 2.0;
			int d;
			for (d=c+1; d<=maxSqsCharge; d++)
			{
				if (sqs_models[d][c].size() <= sqs_size_idx ||
					! sqs_models[d][c][sqs_size_idx])
					continue;

			//	if (c==2 && d == 3)
			//		sqs_sam.print(SQS_var_names);

				comp_prob=sqs_models[d][c][sqs_size_idx]->p_y_given_x(0,sqs_sam);

				if (comp_prob<minComparisonProbs[d])
					minComparisonProbs[d]=comp_prob;

				const float one_minus_prob = 1.0-comp_prob;
				if (one_minus_prob<minComparisonProbs[c])
					minComparisonProbs[c]=one_minus_prob;

			//	cout << "COMP:  " << c << ": " << minComparisonProbs[c] << "\t" << d << ": " << minComparisonProbs[d] << endl;
			}
		}
		
		for (c=1; c<=maximalChargeWithModels_; c++)
			if (minComparisonProbs[c]>1)
				minComparisonProbs[c]=0;
	}
	else // give the max prob charge to the input charge, the rest get 0
	{
		const int charge = pl.getHeader()->getCharge();
		if (charge<=0)
		{
			cout << "Error: no SQS model was read, so charge must be supplied in the spectrum!" << endl;
			exit(1);
		}

		int c;
		for (c=1; c<=maximalChargeWithModels_; c++)
			if (c != charge)
				minComparisonProbs[c]=0;
	}

	
	float max_prob=0;
	int charge;

	const int max_pmc_charge = (pmc_rank_models.size() < maximalChargeWithModels_ ? 
		pmc_rank_models.size() : maximalChargeWithModels_);

	for (charge=1; charge<=max_pmc_charge; charge++)
	{
		const float prob = minComparisonProbs[charge];
		if (prob>max_prob)
			max_prob=prob;

		results[charge].min_comp_prob = prob;

		if (prob>0.02)
		{
		//	find_best_mz_values_from_rank_model(bs,charge,config->get_pm_tolerance(),results[charge]);
			computeBestMzValuesForCharge(pl, charge, config->get_pm_tolerance(), results[charge]);

		//	cout << "RES: " << charge << "\t" << results[charge].score1 << "\t" << results[charge].mz1 << endl;
		//	cout << "RES: " << charge << "\t" << results[charge].score2 << "\t" << results[charge].mz2 << endl;
		}
		else // only give one guess, the second one will come from another charge
		{
			results[charge].mz1= pl.getHeader()->getMOverZ();
			results[charge].score1=0;
			results[charge].mz2=NEG_INF;
			results[charge].score2=NEG_INF;
		}
	}
	

	return max_prob;*/
	return 0.0;
}


bool PMCSQS_Scorer::initializeForCurrentSpectrum(const Config* config, const PeakList& pl)
{
	config_ = config;

	pl.calculateIsotopicLevels(config->getTolerance(), currentSpectrumIsotopeLevels_);
	if (! pl.selectStrongPeakIndexesForPmcsqs(currentSpectrumIsotopeLevels_, currentSpectrumStrongPeakFlags_))
		return false;

	pl.markAllPossibleIsotopicPeaks(config->getTolerance(),currentSpectrumStrictIsotopeFlags_);

	currentSpectrumTotalIntensity_=0;
	currentSpectrumStrongIntensity_=0;
	currentSpectrumNumStrongPeaks_=0;
	const Peak* peaks  = pl.getPeaks();
	const int numPeaks = pl.getNumPeaks();
	for (int i=0; i<numPeaks; i++)
	{
		currentSpectrumTotalIntensity_ += peaks[i].intensity;
		
		if (currentSpectrumIsotopeLevels_[i]>0)
			continue;

		if (currentSpectrumStrongPeakFlags_[i])
		{
			currentSpectrumStrongIntensity_ += peaks[i].intensity;
			currentSpectrumNumStrongPeaks_++;
		}
	}


	if (currentSpectrumPmcTables_.size() <= maximalChargeWithModels_+1)
	{
		currentSpectrumPmcTables_.resize(maximalChargeWithModels_+1);
		currentSpectrumBackgroundStats_.resize(maximalChargeWithModels_+1);
		currentSpectrumMaximalValues_.resize(maximalChargeWithModels_+1);
	}

	return true;
}


void fillPmcRankStatistics(int charge,
                           const mass_t singleChargePairOffset, // the sum of b+y or c+z
                           mass_t minusRange, 
                           mass_t plusRange,
                           mass_t increment,
                           const Config* config,
                           const PeakList& pl,
                           const vector<bool>&   strongPeakIndicators,
                           const vector<float>&  isotopicLevels,
                           const vector<bool>&   strictIsotopicIndicators,
                           vector<PMCRankStats>& pmcStatisticsVector)
{
        const mass_t tolerance   = config->getTolerance()*0.55;
        const int numBinsPerDa   = static_cast<int>(1.0/increment);
        const int numMinusBins   = static_cast<int>((-minusRange) * numBinsPerDa) + 1;
        const int numPlusBins    = static_cast<int>(plusRange * numBinsPerDa) + 1;
        const mass_t oneOverCharge = 1.0/static_cast<mass_t>(charge);

        const int totalNumBins = numMinusBins + numPlusBins;

        if (pmcStatisticsVector.size() != totalNumBins)
                pmcStatisticsVector.resize(totalNumBins);

        int i;
        for (i=0; i<numMinusBins; i++)
        {       
                const mass_t delta = static_cast<mass_t>(i - numMinusBins)*increment;
                calculatePmcRankStatisticsForMass(pl, singleChargePairOffset + delta, tolerance, 
                        isotopicLevels, strictIsotopicIndicators, strongPeakIndicators, pmcStatisticsVector[i]);
                pmcStatisticsVector[i].m_over_z = (singleChargePairOffset + delta + (charge-2)*MASS_PROTON)/
                                                                                  static_cast<mass_t>(charge);
        }

        const int startBinIndex = i;
        for (i=0; i<numPlusBins; i++)
        {
                const mass_t delta = static_cast<mass_t>(i)*increment;

                calculatePmcRankStatisticsForMass(pl, singleChargePairOffset + delta, tolerance, 
                        isotopicLevels, strictIsotopicIndicators, strongPeakIndicators, pmcStatisticsVector[startBinIndex+i]);

                pmcStatisticsVector[startBinIndex+i].m_over_z = (singleChargePairOffset+delta+(charge-2)*MASS_PROTON)/
                                                                                                                  static_cast<mass_t>(charge);
        }
}


void PMCSQS_Scorer::calculatePmcValuesForSqs(const PeakList& pl, mass_t binIncrement)
{
	if (fragmentPairSumOffset < -999.0)
	{
		setFragmentPairSumOffset(MASS_PROTON); // b+y - PM+19
		setBinMassIncrement(0.1);
	}

	const mass_t mOverZ = pl.getHeader()->getMOverZ();
	int charge;
	for (charge=1; charge<=maximalChargeWithModels_; charge++)
	{
		const mass_t thisChargePmWith19 = mOverZ* charge - MASS_PROTON*(charge - 1);
		const mass_t frag_pair_sum = thisChargePmWith19 + fragmentPairSumOffset;

		float bias = 0;
		if (pmcMzBiases_.size()>0 && config_->get_pm_tolerance()>0.075)
		{
			const size_t size_idx = calcRankModelSizeIndex(charge,thisChargePmWith19);
			if (pmcMzBiases_[charge].size()>size_idx)
				bias = pmcMzBiases_[charge][size_idx];
		}

		int bin_range =3*charge;
		if (charge>=2)
			bin_range = 2*charge;
		if (charge>=3)
			bin_range = static_cast<int>(1.5*charge);

		// NP3 GOT BUG wiht currentSpectrumPmcTables_[charge]
		fillPmcRankStatistics(charge,
							  frag_pair_sum + bias,
							  -bin_range-1, 
							  bin_range-1,
							  binIncrement,
							  config_,
							  pl,
							  currentSpectrumStrongPeakFlags_,
							  currentSpectrumIsotopeLevels_,
							  currentSpectrumStrictIsotopeFlags_,
							  currentSpectrumPmcTables_[charge]);
	}
}

void PMCSQS_Scorer::calculateCurrentSpectrumPmcValues(const PeakList& pl, 
													  mass_t binIncrement,
													  int selectedCharge)
{
	if (selectedCharge==0)
		calculatePmcValuesForSqs(pl, binIncrement);

	const mass_t mOverZ = pl.getHeader()->getMOverZ();
	for (int charge=1; charge<=maximalChargeWithModels_; charge++)
	{
		const mass_t thisChargePmWith19 = mOverZ* charge - MASS_PROTON*(charge - 1);
		const mass_t frag_pair_sum = thisChargePmWith19 + fragmentPairSumOffset;
		const size_t size_idx = this->calcRankModelSizeIndex(charge,thisChargePmWith19);

		if (selectedCharge>0 && selectedCharge != charge)
			continue;

		calculateBackgroundStatistics(frag_pair_sum, config_, pl, 
			currentSpectrumStrongPeakFlags_, currentSpectrumIsotopeLevels_, 
			currentSpectrumStrictIsotopeFlags_, currentSpectrumBackgroundStats_[charge]);

		// find maximal values
		int i;
		PMCRankStats& maximal = currentSpectrumMaximalValues_[charge];
		maximal = currentSpectrumPmcTables_[charge][0];

		for (i=1; i<currentSpectrumPmcTables_[charge].size(); i++)
		{
			const PMCRankStats& curr_stats = currentSpectrumPmcTables_[charge][i];

			// mathced pairs look for maximum number

			if (curr_stats.numCharge2FragmentPairs>maximal.numCharge2FragmentPairs)
				maximal.numCharge2FragmentPairs = curr_stats.numCharge2FragmentPairs;

			if (curr_stats.numFragmentPairs>maximal.numFragmentPairs)
				maximal.numFragmentPairs = curr_stats.numFragmentPairs;

			if (curr_stats.numStrongFragmentPairs>maximal.numStrongFragmentPairs)
				maximal.numStrongFragmentPairs = curr_stats.numStrongFragmentPairs;

			if (curr_stats.numStrongCharge2FragmentPairs > maximal.numStrongCharge2FragmentPairs)
				maximal.numStrongCharge2FragmentPairs = curr_stats.numStrongCharge2FragmentPairs;

			// intensity look for maximal numbers

			if (curr_stats.intensityInFragmentPairs > maximal.intensityInFragmentPairs)
				maximal.intensityInFragmentPairs = curr_stats.intensityInFragmentPairs;
			
			if (curr_stats.intensityInStrongPairs> maximal.intensityInStrongPairs)
				maximal.intensityInStrongPairs = curr_stats.intensityInStrongPairs;

			if (curr_stats.intensityInCharge2Pairs> maximal.intensityInCharge2Pairs)
				maximal.intensityInCharge2Pairs = curr_stats.intensityInCharge2Pairs;

			if (curr_stats.intensityInStrongCharge2Pairs > maximal.intensityInStrongCharge2Pairs)
				maximal.intensityInStrongCharge2Pairs = curr_stats.intensityInStrongCharge2Pairs;
		}

		// find indicators for min tolerance and max pairs
		float tol_pairs =           POS_INF;
		float tol_strong_pairs =    POS_INF;
		float tol_c2_pairs =        POS_INF;
		float tol_c2_strong_pairs = POS_INF;

		int   idx_pairs=0;
		int	  idx_strong_pairs=0;
		int   idx_c2_pairs=0;
		int   idx_c2_strong_pairs=0;

		for (i=0; i<currentSpectrumPmcTables_[charge].size(); i++)
		{
			const PMCRankStats& curr_stats = currentSpectrumPmcTables_[charge][i];

			if (curr_stats.numFragmentPairs == maximal.numFragmentPairs &&
				curr_stats.meanOffsetPairs < tol_pairs)
			{
				idx_pairs=i;
				tol_pairs=curr_stats.meanOffsetPairs;
			}

			if (curr_stats.numStrongFragmentPairs == maximal.numStrongFragmentPairs &&
				curr_stats.meanOffsetCharge2StrongPairs < tol_strong_pairs)
			{
				idx_strong_pairs=i;
				tol_strong_pairs=curr_stats.meanOffsetCharge2StrongPairs;
			}

			if (curr_stats.numCharge2FragmentPairs == maximal.numCharge2FragmentPairs &&
				curr_stats.meanOffsetCharge2Pairs < tol_c2_pairs)
			{
				idx_c2_pairs=i;
				tol_c2_pairs=curr_stats.meanOffsetCharge2Pairs;
			}

			if (curr_stats.numStrongCharge2FragmentPairs == maximal.numStrongCharge2FragmentPairs &&
				curr_stats.meanOffsetCharge2StrongPairs < tol_c2_strong_pairs)
			{
				idx_c2_strong_pairs=i;
				tol_c2_strong_pairs=curr_stats.meanOffsetCharge2StrongPairs;
			}
		}

		currentSpectrumPmcTables_[charge][idx_pairs].ind_pairs_with_min_tol=true;
		currentSpectrumPmcTables_[charge][idx_strong_pairs].ind_strong_pairs_with_min_tol=true;
		currentSpectrumPmcTables_[charge][idx_c2_pairs].ind_c2_pairs_with_min_tol=true;
		currentSpectrumPmcTables_[charge][idx_c2_strong_pairs].ind_c2_strong_pairs_with_min_tol=true;

	}
}



PMCRankStats::PMCRankStats(const PMCRankStats & Source)
{
	*this = Source;
//	offsetPairsOderedByIntensity		 = Source.offsetPairsOderedByIntensity;
//	strong_offsetPairsOderedByIntensity = Source.strong_offsetPairsOderedByIntensity;
//  charge2offsetPairsOderedByIntensity	 =  Source.charge2offsetPairsOderedByIntensity;
}


void PMCSQS_Scorer::set_sqs_mass_thresholds()
{
	sqsMassThresholds_.clear();
	sqsMassThresholds_.push_back(800);
	sqsMassThresholds_.push_back(1200);
}










/***************************************************************************
Gives detailed m/z and prob values for the different charges.
Returns the highest sqs probability found for a spectrum.
****************************************************************************
float PMCSQS_Scorer::get_pmcsqs_results_for_spectrum(const Config *config, 
													 const BasicSpectrum& bs,
													 vector<PmcSqsChargeRes>& res)
{
	static ME_Regression_Sample sqs_sam;

	if (! init_for_current_spec(config,bs))
		return 0; // corrupt file

	calculate_curr_spec_pmc_values(bs,bin_increment);

	const int sqs_size_idx = this->getSqsSizeIndex(bs.ssf->m_over_z);

	if (res.size()<=maximalChargeWithModels_)
		res.resize(maximalChargeWithModels_+1);

	
	vector<float> minComparisonProbs;		// the minimal probability when used in a comparison model
	minComparisonProbs.resize(maximalChargeWithModels_+1,2.0);
	minComparisonProbs[0]=0;

	if (indInitializedSqs_)
	{
		const int maxSqsCharge = ( sqs_models.size()<maximalChargeWithModels_ ?
			sqs_models.size() : maximalChargeWithModels_);

	//	fill_fval_vector_with_SQS(bs, sqs_sam);
	
		int charge;
		for (charge=1; charge<=maxSqsCharge; charge++)
		{
			if (sqs_models[charge].size()<1 ||
				sqs_models[charge][0].size() <= sqs_size_idx)
				continue;
			
			float prob = sqs_models[charge][0][sqs_size_idx]->p_y_given_x(0,sqs_sam);

			// correct prob
			prob += sqs_correction_factors[charge][sqs_size_idx];
			prob *= sqs_mult_factors[charge][sqs_size_idx];
			if (prob<0.0)
				prob=0.0;

			res[charge].sqs_prob=prob;
		}
		
		int c;
		for (c=1; c<maxSqsCharge; c++)
		{
			if (sqs_models[c].size()<=0)
				continue;

			float comp_prob = 2.0;
			int d;
			for (d=c+1; d<=maxSqsCharge; d++)
			{
				if (sqs_models[d][c].size() <= sqs_size_idx ||
					! sqs_models[d][c][sqs_size_idx])
					continue;

				comp_prob=sqs_models[d][c][sqs_size_idx]->p_y_given_x(0,sqs_sam);

				if (comp_prob<minComparisonProbs[d])
					minComparisonProbs[d]=comp_prob;

				const float one_minus_prob = 1.0-comp_prob;
				if (one_minus_prob<minComparisonProbs[c])
					minComparisonProbs[c]=one_minus_prob;
			}
		}
		
		for (c=1; c<=maximalChargeWithModels_; c++)
			if (minComparisonProbs[c]>1)
				minComparisonProbs[c]=0;
	}
	else // give the max prob charge to the input charge, the rest get 0
	{
		if (bs.ssf->charge<=0)
		{
			cout << "Error: no SQS model was read, so charge must be supplied in the spectrum!" << endl;
			exit(1);
		}

		int c;
		for (c=1; c<=maximalChargeWithModels_; c++)
			if (c != bs.ssf->charge)
				minComparisonProbs[c]=0;
	}

	
	float max_prob=0;
	int charge;

	const int max_pmc_charge = (pmc_rank_models.size() < maximalChargeWithModels_ ? 
		pmc_rank_models.size() : maximalChargeWithModels_);

	for (charge=1; charge<=max_pmc_charge; charge++)
	{
		const float prob = minComparisonProbs[charge];
		if (prob>max_prob)
			max_prob=prob;

		res[charge].min_comp_prob = prob;

		if (prob>0.02)
		{
			find_best_mz_values_from_rank_model(bs,charge,config->get_pm_tolerance(),res[charge]);
		}
		else // only give one guess, the second one will come from another charge
		{
			res[charge].mz1=bs.ssf->m_over_z;
			res[charge].score1=0;
			res[charge].mz2=NEG_INF;
			res[charge].score2=NEG_INF;
		}
	}
	

	return max_prob;
	return 0.0;
}
*/




/********************************************************************
Computes the best sqs value for the spectrum. It is faster than the charge
m/z function since it only looks at a few mass positions. 
*********************************************************************
float PMCSQS_Scorer::get_sqs_for_spectrum(const Config *config, 
										  const BasicSpectrum& bs, 
										  int *max_charge,
										  bool verbose)
{
	static ME_Regression_Sample sqs_sam;
	static vector<QCPeak> sqs_peaks;
	static int num_sqs_peaks = 0;
	BasicSpectrum sqs_bs;

	if (verbose)
	{
		cout << "ORG spectrum: " << bs.num_peaks << endl;
//		bs.print_peaks();
		cout << endl;
	}

	if (num_sqs_peaks< 2000 || num_sqs_peaks<bs.num_peaks)
	{
		num_sqs_peaks = (int)(bs.num_peaks * 1.5);
		if (num_sqs_peaks<2000)
			num_sqs_peaks = 2000;
		
		if (num_sqs_peaks>100000)
		{
			cout << "Error: too many peak in spectrum: " << bs.num_peaks << endl;
			exit(1);
		}
		sqs_peaks.resize(num_sqs_peaks);
	}

	sqs_bs.peaks = &sqs_peaks[0];
	sqs_bs.ssf = bs.ssf;

	int best_charge=0;

	// filter peaks to acheive the required peak density for the models

	int new_num_peaks =0;
	create_filtered_peak_list_for_sqs(bs.peaks,bs.num_peaks,sqs_bs.peaks,new_num_peaks);
	sqs_bs.num_peaks= new_num_peaks;

	if (verbose)
	{
		cout << "ATER FILTERING: " << new_num_peaks << endl;
//		sqs_bs.print_peaks();
		cout << endl;
	}

	if (! init_for_current_spec(config,sqs_bs))
		return 0; // corrupt spectrum

	calculate_curr_spec_pmc_values(sqs_bs,bin_increment*3);

	const int num_charges = pmc_rank_models.size();
	const int sqs_size_idx = this->getSqsSizeIndex(sqs_bs.ssf->m_over_z);

	float max_prob=-1;
	
	//fill_fval_vector_with_SQS(sqs_bs, sqs_sam);
	
	int charge;
	for (charge=1; charge<sqs_models.size(); charge++)
	{
		if (sqs_models[charge].size()<1 ||
			sqs_models[charge][0].size()<=sqs_size_idx ||
			! sqs_models[charge][0][sqs_size_idx] ||
			! sqs_models[charge][0][sqs_size_idx]->get_has_weights() )
			continue;
		
		float prob = sqs_models[charge][0][sqs_size_idx]->p_y_given_x(0,sqs_sam);

		// correct prob
		prob += sqs_correction_factors[charge][sqs_size_idx];
		prob *= sqs_mult_factors[charge][sqs_size_idx];
		if (prob<0.0)
			prob=0.0;

		if (prob>max_prob)
		{
			max_prob=prob;
			best_charge=charge;
		}
	}

	if (max_charge)
		*max_charge = best_charge;

	return max_prob;
	return 0.0;
}*/



float PMCSQS_Scorer::computeBestMzAndChargeForSpectrum(const Config *config, const PeakList& pl,
							mass_t* mz1, int* charge1, float *prob1,
							mass_t* mz2, int* charge2, float *prob2,
							vector<PmcSqsChargeRes>* all_res)
{
	static vector<PmcSqsChargeRes> res;

	if (! this->indInitializedPmc_)
	{
		cout << "Error: no PMC model was read!" << endl;
		if (pl.getHeader()->getCharge() <=0)
			cout << "Spectrum was read with charge <=0, so charge selection is needed!" << endl;
		exit(1);
	}

//	get_pmcsqs_results_for_spectrum(config,bs,res);
	computePmcsqsResultsForSpectrum(config, pl, res);

	float best_prob=-1;
	int best_charge=0;

	int charge;
	for (charge =1; charge<=maximalChargeWithModels_; charge++)
		if (res[charge].min_comp_prob>best_prob)
		{
			best_charge = charge;
			best_prob = res[charge].min_comp_prob;
		}

	const PmcSqsChargeRes& cr = res[best_charge];
	float second_best_prob=1.0 - (cr.score1-cr.score2)/(fabs(cr.score1)+fabs(cr.score2)+2.0);
	int   second_best_charge = best_charge;
	
	for (charge =1; charge<=maximalChargeWithModels_; charge++)
		if (charge != best_charge && res[charge].min_comp_prob>second_best_prob)
		{
			second_best_charge = charge;
			second_best_prob = res[charge].min_comp_prob;
		}

	*mz1 = (mass_t)res[best_charge].mz1;
	*charge1 = best_charge;
	*prob1 = best_prob;

	if (mz2 && charge2)
	{
		*charge2 = second_best_charge;
		*prob2   = second_best_prob;
		if (second_best_charge == best_charge)
		{
			*mz2 = (mass_t)res[best_charge].mz2;
		}
		else
			*mz2 = (mass_t)res[second_best_charge].mz1;
	}

	if (all_res)
		*all_res = res;

	return best_prob;
}

/******************************************************************************
Selects the the two best values of charges 1,2,3
returns the max prob found
*******************************************************************************
float PMCSQS_Scorer::get_best_mz_charge(const Config *config, const BasicSpectrum& bs, 
						   mass_t* mz1, int* charge1, float *prob1,
						   mass_t* mz2, int* charge2, float *prob2,
						   vector<PmcSqsChargeRes>* all_res)
{
	static vector<PmcSqsChargeRes> res;

	if (! this->indInitializedPmc_)
	{
		cout << "Error: no PMC model was read!" << endl;
		if (bs.ssf->charge <=0)
			cout << "Spectrum was read with charge <=0, so charge selection is needed!" << endl;
		exit(1);
	}

//	get_pmcsqs_results_for_spectrum(config,bs,res);

	float best_prob=-1;
	int best_charge=0;

	int charge;
	for (charge =1; charge<=maximalChargeWithModels_; charge++)
		if (res[charge].min_comp_prob>best_prob)
		{
			best_charge = charge;
			best_prob = res[charge].min_comp_prob;
		}

	const PmcSqsChargeRes& cr = res[best_charge];
	float second_best_prob=1.0 - (cr.score1-cr.score2)/(fabs(cr.score1)+fabs(cr.score2)+2.0);
	int   second_best_charge = best_charge;
	
	for (charge =1; charge<=maximalChargeWithModels_; charge++)
		if (charge != best_charge && res[charge].min_comp_prob>second_best_prob)
		{
			second_best_charge = charge;
			second_best_prob = res[charge].min_comp_prob;
		}

	*mz1 = (mass_t)res[best_charge].mz1;
	*charge1 = best_charge;
	*prob1 = best_prob;

	if (mz2 && charge2)
	{
		*charge2 = second_best_charge;
		*prob2   = second_best_prob;
		if (second_best_charge == best_charge)
		{
			*mz2 = (mass_t)res[best_charge].mz2;
		}
		else
			*mz2 = (mass_t)res[second_best_charge].mz1;
	}

	if (all_res)
		*all_res = res;

	return best_prob;
}*/




/*
void PMCSQS_Scorer::benchmark_pm_selection(Config *config, FileManager& fm, mass_t pm_val_tol)
{
	const vector< vector< mass_t > >& threshes = config->get_size_thresholds();
	int c=1;

	for (c=1; c<threshes.size(); c++)
	{
		int s;
		for (s=0; s<threshes[c].size(); s++)
		{
			mass_t min_mz = (s>0 ? threshes[c][s-1]/c : 0);
			mass_t max_mz = threshes[c][s]/c;

			FileSet fs;
			fs.select_files_in_mz_range(fm,min_mz,max_mz,c);

			if (fs.get_total_spectra()<200)
				continue;

			cout << "CHARGE " << c <<" size " << s << "  (" << fs.get_total_spectra() << " spectra)" << endl;

			const vector<SingleSpectrumFile *>& all_ssfs = fs.get_ssf_pointers();
			BasicSpecReader bsr;
			vector<QCPeak> peaks;
			peaks.resize(10000);

			vector<int> correct_counts;
			correct_counts.resize(8,0);

			int num_correct=0;
			int num_wrong_charge=0;
			int num_diff_charge=0;

			int i;
			for (i=0; i<all_ssfs.size(); i++)
			{
				SingleSpectrumFile *ssf=all_ssfs[i];
				BasicSpectrum bs;

				bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,&peaks[0]);
				bs.peaks = &peaks[0];
				bs.ssf = ssf;

				const mass_t true_mass = ssf->peptide.get_mass_with_19();

				vector<mass_t> pms_with_19;
				vector<int>    charges;
				select_pms_and_charges(config,bs,pms_with_19,charges);

				if (charges[0] != c)
					num_wrong_charge++;

				bool got_diff_charge=false;
				int j;
				for (j=1; j<charges.size(); j++)
					if (charges[j] != charges[0])
						got_diff_charge=true;

				if (got_diff_charge)
					num_diff_charge++;

				for (j=0; j<pms_with_19.size(); j++)
					if (fabs(pms_with_19[j]-true_mass)<pm_val_tol)
						break;
				
				if (j==pms_with_19.size())
					continue;
				

				num_correct++;

				if (got_diff_charge && j == pms_with_19.size()-1)
				{
					correct_counts[7]++;
				}
				else
					correct_counts[j]++;
			}

			double num_total = (double)all_ssfs.size();
			cout << "Had correct       " << fixed << setprecision(4) << num_correct/num_total << endl;
			cout << "First correct     " << correct_counts[0]/num_total << endl;
			cout << "Second correct    " << correct_counts[1]/num_total << endl;
			cout << "Off-1  correct    " << correct_counts[2]/num_total << endl;
			cout << "Off+1  correct    " << correct_counts[3]/num_total << endl;
			cout << "Off-1  2nd        " << correct_counts[4]/num_total << endl;
			cout << "Off+1  2nd        " << correct_counts[5]/num_total << endl;
			cout << "Diff Ch correct   " << correct_counts[7]/num_total << endl;
			cout << "With wrong charge " << num_wrong_charge/num_total << endl;
			cout << "With diff charge  " << num_diff_charge/num_total << endl << endl;
		}
	}
}

*/





// for each peak,aa holds the idx of the previous peak if they have a mass
// diff of that aa (within tolerance)
// entry 0 in each column holds an indicator if peak is in aa diff
// entry 1 in each column holds an indicator if peak has a 
void PMCSQS_Scorer::fillSqsDynamicProgrammingTable(const PeakList& pl, 
												   vector<DPColumn>& dp, 
												   int fargmentCharge ) const
{
	const Peak* const peaks = pl.getPeaks();
	const size_t numPeaks   = pl.getNumPeaks();
	const vector<mass_t>& aa2mass = config_->get_aa2mass();
	const mass_t tagTolerance    = (config_->getTolerance()<0.15 ? 
									config_->getTolerance() : config_->getTolerance()*0.33);

	const mass_t multVal = (1.0 / static_cast<mass_t>(fargmentCharge));
	dp.resize(numPeaks);
	int i;
	for (i=0; i<numPeaks; i++)
	{
		dp[i].pointers[0]=0;
		dp[i].fPointers[0]=0;
		int j;
		for (j=1; j<=Val; j++)
		{
			dp[i].pointers[j]=-1;
			dp[i].fPointers[j]=-1;
		}
	}

	int aa;
	for (aa=Ala; aa<=Val; aa++)
	{
		if (aa==Ile || aa==Xle)
			continue;

		const mass_t aaMass = multVal * aa2mass[aa];
		const mass_t minimalOffset = aaMass-tagTolerance;
		const mass_t maximalOffset = aaMass+tagTolerance;

		size_t trailIndex=0;
		size_t leadIndex=1;

		while (leadIndex<numPeaks)
		{
			if (currentSpectrumIsotopeLevels_[leadIndex]>0)
			{
				leadIndex++;
				continue;
			}

			while (peaks[leadIndex].mass-peaks[trailIndex].mass>maximalOffset)
				trailIndex++;

			if (currentSpectrumIsotopeLevels_[trailIndex]==0 && 
				peaks[leadIndex].mass-peaks[trailIndex].mass>minimalOffset)
			{
				dp[leadIndex].pointers[aa]=trailIndex;
				dp[leadIndex].pointers[0]=1;
				dp[trailIndex].pointers[0]=1;
				dp[trailIndex].fPointers[aa]=leadIndex;
			}
			leadIndex++;
		}
	}
}





// used in fillSqsSample for sorting peaks according to increasing intensity
struct IntenMassPair
{
	bool operator< (const IntenMassPair& rhs) const
	{
		return (inten>rhs.inten);
    }
    intensity_t inten;
    mass_t      mass;
};



/**************************************************************************
Fills a sample with the SQS features (listed in sqsFeatures.txt)
***************************************************************************/
void PMCSQS_Scorer::fillSqsSample(const PeakList& pl, MlSample& sample) const
{
	const mass_t tolerance = config_->getTolerance();
	const mass_t fragmentTolerance = (tolerance<0.15 ? tolerance : tolerance * 0.5);
	const Peak* const peaks  = pl.getPeaks();
	const int  numPeaks		 = pl.getNumPeaks();
	const mass_t mOverZ		 = pl.getHeader()->getMOverZ();
	const int numStrongPeaks = currentSpectrumNumStrongPeaks_;
	const float totalIntensity = currentSpectrumTotalIntensity_;
	const float oneOverNumberOfPeaks = 1.0 / (numPeaks+1);
	const float oneOverNumberOfStrongPeaks = 1.0/ (numStrongPeaks + 1);
	const float oneOverTotalIntensity = 1.0 / (totalIntensity + 1.0);
	const mass_t maximalPeakMass = (numPeaks <5 ? POS_INF : peaks[numPeaks-1].mass);


	sample.clearPairs();

	sample.addPair(0,1.0); // SQS_CONST
	sample.addPair(1, static_cast<value_t>(numPeaks)/maximalPeakMass); // SQS_PEAK_DENSITY

	float grassLevelIntensity=0;
	
	// calculate grass level peaks and add intensity features
	if (1)
	{
		int i;
		vector<float> peakIntensities;
		peakIntensities.resize(numPeaks);
		for (i=0; i<numPeaks; i++)
			peakIntensities[i]=peaks[i].intensity;

		sort(peakIntensities.begin(),peakIntensities.end());

		size_t indexForGrass = numPeaks/3;
		grassLevelIntensity = peakIntensities[indexForGrass];

		float intensityUpto1G=0.0;
		size_t index=0;
		while (index<numPeaks && peakIntensities[index]<grassLevelIntensity)
			intensityUpto1G+=peakIntensities[index++];

		float intensityUpto2G=0.0;
		float intensityOf2G = 2.0 * grassLevelIntensity;
		size_t index2G = indexForGrass;
		while (index2G<numPeaks && peakIntensities[index2G]<intensityOf2G)
			intensityUpto2G+=peakIntensities[index2G++];
		
		float intensityUpto5G=0.0;
		float intensityOf5G = 5.0 * grassLevelIntensity;
		size_t index5G = index2G;
		while (index5G<numPeaks && peakIntensities[index5G]<intensityOf5G)
			intensityUpto5G+=peakIntensities[index5G++];
		
		float intensityUpto10G=0.0;
		float intensityOf10G = 10.0 * grassLevelIntensity;
		size_t index10G = index5G;
		while (index10G<numPeaks && peakIntensities[index10G]<intensityOf10G)
			intensityUpto10G+=peakIntensities[index10G++];

		if (indexForGrass>0)
			sample.addPair(2, static_cast<value_t>(indexForGrass)* oneOverNumberOfPeaks);  // SQS_PROP_GRASS

		if (index2G-indexForGrass>0)
			sample.addPair(3, static_cast<value_t>(index2G-indexForGrass)* oneOverNumberOfPeaks); // SQS_PROP_PEAKS_1GTO2G

		if (index5G-index2G>0)
			sample.addPair(4, static_cast<value_t>(index5G-index2G)* oneOverNumberOfPeaks); // SQS_PROP_PEAKS_2GTO5G

		if (index10G-index5G>0)
			sample.addPair(5, static_cast<value_t>(index10G-index5G)* oneOverNumberOfPeaks); // SQS_PROP_PEAKS_5GTO10G

		if (numPeaks-index10G>0)
			sample.addPair(6, static_cast<value_t>(numPeaks-index10G)* oneOverNumberOfPeaks); // SQS_PROP_PEAKS_MORE_10G

		if (intensityUpto1G>0.0)
			sample.addPair(7, intensityUpto1G * oneOverTotalIntensity); // SQS_PROP_INTENSITY_GRASS

		if (intensityUpto2G>0.0)
			sample.addPair(8, intensityUpto2G * oneOverTotalIntensity); // SQS_PROP_INTENSITY_2G

		if (intensityUpto5G>0.0)
			sample.addPair(9, intensityUpto5G * oneOverTotalIntensity); // SQS_PROP_INTENSITY_5G

		if (intensityUpto10G>0.0)
			sample.addPair(10, intensityUpto10G * oneOverTotalIntensity); // SQS_PROP_INTENSITY_10G

		const float intenDiff = totalIntensity-intensityUpto1G-intensityUpto2G-intensityUpto5G-intensityUpto10G;
		if (intenDiff>0.0)
		sample.addPair(11, intenDiff * oneOverTotalIntensity); // SQS_PROP_INTENSITY_MORE_10G
	}

	// isotope features
	if (1)
	{
		size_t i;
		size_t numPeaksWithIso=0;
		size_t numStrongPeaksWithIso=0;
		for (i=1; i<numPeaks; i++)
			if (currentSpectrumIsotopeLevels_[i]>0)
			{
				numPeaksWithIso++;
				if (currentSpectrumIsotopeLevels_[i-1]==0.0 && currentSpectrumStrongPeakFlags_[i-1])
					numStrongPeaksWithIso++;
			}

		sample.addPair(12, numPeaksWithIso*oneOverNumberOfPeaks); // SQS_PROP_PEAKS_WITH_ISO
		sample.addPair(13, numStrongPeaksWithIso*oneOverNumberOfStrongPeaks); // PROP_STRONG_PEAKS_WITH_ISO
	}


	// neutral loss features
	if (1)
	{
		int anyLoss[3]={0};
		int anyStrongLoss[3]={0};
		for (size_t charge=1; charge<=2; charge++)
		{
			const mass_t offsets[]={MASS_H2O/charge, MASS_NH3/charge, MASS_CO/charge};
			const size_t numOffsets = sizeof(offsets)/sizeof(mass_t);
			int	  numPairs[numOffsets]={0};
			int	  numStrongPairs[numOffsets]={0};

			for (size_t i=0; i<numOffsets; i++)
			{
				// find how many pairs of peaks and offset exist
				// (also look at how many of them involve a 
				const mass_t minimalOffset = offsets[i] - fragmentTolerance;
				const mass_t maximalOffset = offsets[i] + fragmentTolerance;
				size_t trailIndex=0, leadIndex=1;
				while (leadIndex<numPeaks)
				{
					while (peaks[leadIndex].mass-peaks[trailIndex].mass>maximalOffset)
						trailIndex++;

					if (peaks[leadIndex].mass-peaks[trailIndex].mass>minimalOffset)
					{
						numPairs[i]++;
						if (currentSpectrumStrongPeakFlags_[leadIndex])
							numStrongPairs[i]++;
					}
					leadIndex++;
				}
				anyLoss[charge] += numPairs[i]; // don't care if same peak gets recounted
				anyStrongLoss[charge] += numStrongPairs[i];
			}
			const size_t indexStart = 14+(charge-1)*6;
			if (numPairs[0]>0)
				sample.addPair(indexStart,  numPairs[0]);  // SQS_NUM_C1_PEAKS_H2O_LOSS
			if (numStrongPairs[0]>0)
				sample.addPair(indexStart+1, numStrongPairs[0]); // SQS_NUM_C1_STRONG_H2O_LOSS
			if (numPairs[1]>0)
				sample.addPair(indexStart+2,  numPairs[1]);  // SQS_NUM_C1_PEAKS_NH3_LOSS
			if (numPairs[2]>0)
				sample.addPair(indexStart+3,  numPairs[2]);  // SQS_NUM_C1_PEAKS_CO_LOSS
			if (anyLoss[charge]>0)
				sample.addPair(indexStart + 4, anyLoss[charge]); // ANY_LOSS
			if (anyStrongLoss[charge]>0)
				sample.addPair(indexStart + 5, anyStrongLoss[charge]); // ANY STRONG LOSS 
		}

		// sum and diff features
		sample.addPair(26, anyLoss[1] + anyLoss[2]);			 // SQS_SUM_C1_C2_ANY_LOSS
		sample.addPair(27, anyStrongLoss[1] + anyStrongLoss[2]); // SQS_SUM_C1_C2_STRONG_LOSS

		const float maxAnyLoss = (anyLoss[1]>anyLoss[2] ? 
									static_cast<float>(anyLoss[1]) : static_cast<float>(anyLoss[2]));
		if (maxAnyLoss>0.0)
			sample.addPair(28, (anyLoss[1] - anyLoss[2])/ maxAnyLoss); // SQS_DIFF_C1_C2_ANY_LOSS

		const float maxAnyStrongLoss = (anyStrongLoss[1]>anyStrongLoss[2] ? 
									static_cast<float>(anyStrongLoss[1]) : static_cast<float>(anyStrongLoss[2]));
		if (maxAnyStrongLoss>0.0)
			sample.addPair(29, (anyStrongLoss[1] - anyStrongLoss[2])/ maxAnyStrongLoss); // SQS_DIFF_C1_C2_STRONG_ANY_LOSS
	}
		
	//	double/single complement features
	if (1)
	{
		const mass_t halfMaxMass = maximalPeakMass*0.5;
		const mass_t smallerTolerance = tolerance * 0.66;
		
		size_t numPairs=0;
		size_t numStrongPairs=0;
		float pairIntensity=0;
		size_t doublePeakIndex=0;
		for (size_t i=0; i<numPeaks; i++)
		{
			const mass_t peakMass = peaks[i].mass;
			if (peakMass>halfMaxMass)
				break;

			const mass_t doublePeakMass = peakMass * 2 - MASS_PROTON;
			const mass_t minimalMass = doublePeakMass - smallerTolerance;
			const mass_t maximalMass = doublePeakMass + smallerTolerance;
			while (doublePeakIndex<numPeaks && peaks[doublePeakIndex].mass<minimalMass)
				doublePeakIndex++;

			if (doublePeakIndex==numPeaks)
				break;
			
			if (peaks[doublePeakIndex].mass<maximalMass)
			{
				numPairs++;
				if (currentSpectrumStrongPeakFlags_[i] || currentSpectrumStrongPeakFlags_[doublePeakIndex])
					numStrongPairs++;
				pairIntensity+= peaks[i].intensity + peaks[doublePeakIndex].intensity;
			}
		}

		sample.addPair(30, numPairs);	//	SQS_PROP_PEAKS_WITH_C1C2_PAIRS
		sample.addPair(31, numStrongPairs); // SQS_PROP_STRONG_WITH_C1C2_PAIRS
		sample.addPair(32, pairIntensity*oneOverTotalIntensity); // SQS_PROP_INTEN_IN_C1C2_PAIRS
	}

	// tag features
	if (1)
	{
		vector<DPColumn> dp;
		vector<float> tmpVals;
		int charge;
		for (charge=1; charge<=2; charge++)
		{
			fillSqsDynamicProgrammingTable(pl, dp, charge);

			vector<int> forwardTagLengths(numPeaks,0);
			vector<int> backwardsTagLengths(numPeaks,0);

			for (int i=1; i<numPeaks; i++)
			{
				if (dp[i].pointers[0]<1)
					continue;
				int maxTagLength = -1;
				for (int aa=Ala; aa<=Val; aa++)
				{
					const int& prev = dp[i].pointers[aa];
					if (prev>=0 && backwardsTagLengths[prev]>maxTagLength)
						maxTagLength = backwardsTagLengths[prev];
				}
				backwardsTagLengths[i]=maxTagLength+1;
			}

			for (int i=numPeaks-1; i>0 ; i--)
			{
				if (dp[i].pointers[0]<1)
					continue;

				int maxTagLength = -1;
				for (int aa=Ala; aa<=Val; aa++)
				{
					const int& prev = dp[i].fPointers[aa];
					if (prev>=0 && forwardTagLengths[prev]>maxTagLength)
						maxTagLength = forwardTagLengths[prev];
				}
				forwardTagLengths[i]=maxTagLength+1;
			}

			vector<int> strongCountsPerLength(6,0);
			int maxTagLength=0;
			for (int i=0; i<numPeaks; i++)
			{
				if (! currentSpectrumStrongPeakFlags_[i])
					continue;
				int len = backwardsTagLengths[i]+forwardTagLengths[i];
				if (len<0)
					continue;

				if (len>maxTagLength)
					maxTagLength=len;

				if (len>=5)
					len=5;
				strongCountsPerLength[len]++;
			}

			const size_t startIndex = (charge == 1 ? 33 : 42);

			if (maxTagLength<=2)
			{
				sample.addPair(startIndex, 1.0);
			}
			else if (maxTagLength<=4)
			{
				sample.addPair(startIndex+1, 1.0);
			}
			else if (maxTagLength<=6)
			{
				sample.addPair(startIndex+2, 1.0);
			}
			else
				sample.addPair(startIndex+3, 1.0);

			if (strongCountsPerLength[0]>0)
				sample.addPair(startIndex+4, strongCountsPerLength[0]);

			if (strongCountsPerLength[3]>0)
			{
				sample.addPair(startIndex+5, strongCountsPerLength[3]);
				sample.addPair(startIndex+6, strongCountsPerLength[3]-strongCountsPerLength[0]);
			}

			if (strongCountsPerLength[5]>0)
			{
				sample.addPair(startIndex+7, strongCountsPerLength[5]);
				sample.addPair(startIndex+8, strongCountsPerLength[5]-strongCountsPerLength[0]);
			}
		}
	}


	// regional density features
	if (1)
	{
		mass_t t1=mOverZ*0.6666;
		mass_t t2=mOverZ*1.3333;
		mass_t h1=mOverZ;

		float intensityT1=0,intensityT2=0, intensityH1=0;
		size_t numPeaksT1=0, numPeaksT2=0, numPeaksH1=0;

		size_t i;
		for (i=0; i<numPeaks && peaks[i].mass<t1; i++)
			intensityT1+=peaks[i].intensity;
		
		numPeaksT1=i;

		for ( ; i<numPeaks && peaks[i].mass<h1; i++)
			intensityH1+=peaks[i].intensity;

		numPeaksH1=i;

		intensityT2=intensityH1;

		for ( ; i<numPeaks && peaks[i].mass<t2; i++)
			intensityT2+=peaks[i].intensity;

		numPeaksT2=i-numPeaksT1;

		const size_t numPeaksT3 = numPeaks - numPeaksT2 - numPeaksT1;
		const float intensityT3 = totalIntensity - intensityT2 - intensityT1;

		sample.addPair(51, static_cast<value_t>(numPeaksT1));
		sample.addPair(52, static_cast<value_t>(numPeaksT2));
		sample.addPair(53, static_cast<value_t>(numPeaksT3));

		sample.addPair(54,	intensityT1 * oneOverTotalIntensity);
		sample.addPair(55,	intensityT2 * oneOverTotalIntensity);
		sample.addPair(56,	intensityT3 * oneOverTotalIntensity);

		const size_t numPeaksH2 = numPeaks - numPeaksH1;
		const float intensityH2 = totalIntensity - intensityH1;

		sample.addPair(57, static_cast<value_t>(numPeaksH1));
		sample.addPair(58, static_cast<value_t>(numPeaksH2));
		sample.addPair(59,	intensityH1 * oneOverTotalIntensity);
		sample.addPair(60,	intensityH2 * oneOverTotalIntensity);
	
		const float inten33=0.333* totalIntensity;
		const float inten50=0.5  * totalIntensity;
		const float inten75=0.75 * totalIntensity;
		const float inten90=0.90 * totalIntensity;

		float cumulativeIntensity=0;
		const mass_t oneOverMz = 1.0 / (mOverZ + 0.1);

		for (i=0; i<numPeaks && cumulativeIntensity<inten33; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==numPeaks)
			i--;

		sample.addPair(61, (peaks[i].mass*oneOverMz)-0.333 ); // SQS_PROP_MZ_RANGE_WITH_33_INTEN

		for ( ; i<numPeaks && cumulativeIntensity<inten50; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==numPeaks)
			i--;

		sample.addPair(62, (peaks[i].mass*oneOverMz)-0.5 ); //  SQS_PROP_MZ_RANGE_WITH_50_INTEN


		for ( ; i<numPeaks && cumulativeIntensity<inten75; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==numPeaks)
			i--;

		sample.addPair(63, (peaks[i].mass*oneOverMz)-0.75 ); //  SQS_PROP_MZ_RANGE_WITH_75_INTEN

		for ( ; i<numPeaks && cumulativeIntensity<inten90; i++)
			cumulativeIntensity+=peaks[i].intensity;

		if (i==numPeaks)
			i--;

		sample.addPair(64, (peaks[i].mass*oneOverMz)-0.9 ); //  SQS_PROP_MZ_RANGE_WITH_90_INTEN	
	}

	// pmc features
	if (1)
	{
		vector< vector<value_t> > sqsFeatures;
		getSqsFeaturesFromPmcTabels(pl, sqsFeatures);

		// find what charges give the best fragment pair counts
		int maxCharge=0;
		for (int charge=1; charge<=maximalChargeWithModels_; charge++)
			if (sqsFeatures[charge][1]>=sqsFeatures[maxCharge][1])
				maxCharge=charge;

		if (maxCharge == 1)
		{
			sample.addPair(65, 1.0); // SQS_IND_CHARGE_MAX_STRONG_PAIRS_IS_1
		}
		else if (maxCharge == 2)
		{
			sample.addPair(66, 1.0); // SQS_IND_CHARGE_MAX_STRONG_PAIRS_IS_2
		}
		else if (maxCharge == 3)
		{
			sample.addPair(67, 1.0); // SQS_IND_CHARGE_MAX_STRONG_PAIRS_IS_3
		}
		else if (maxCharge > 3)
			sample.addPair(68, 1.0); // SQS_IND_CHARGE_MAX_STRONG_PAIRS_IS_4+

		// find what charges give the best fragment pair c2 strong counts
		maxCharge=0;
		for (int charge=1; charge<=maximalChargeWithModels_; charge++)
			if (sqsFeatures[charge][3]>=sqsFeatures[maxCharge][3])
				maxCharge=charge;

		if (maxCharge == 1)
		{
			sample.addPair(69, 1.0); // SQS_IND_CHARGE_MAX_C2_STRONG_PAIRS_IS_1
		}
		else if (maxCharge == 2)
		{
			sample.addPair(70, 1.0); // SQS_IND_CHARGE_MAX_C2_STRONG_PAIRS_IS_2
		}
		else if (maxCharge == 3)
		{
			sample.addPair(71, 1.0); // SQS_IND_CHARGE_MAX_C2_STRONG_PAIRS_IS_3
		}
		else if (maxCharge > 3)
			sample.addPair(72, 1.0); // SQS_IND_CHARGE_MAX_C2_STRONG_PAIRS_IS_4+

		const int maxChargeForLoop = (maximalChargeWithModels_ <=4 ? maximalChargeWithModels_ : 4);
		for (int charge=1; charge<=maxChargeForLoop; charge++)
			sample.addPair(72+charge, sqsFeatures[charge][0]); // SQS_NUM_PAIRS_CHARGE

		for (int charge=1; charge<=maxChargeForLoop; charge++)
			sample.addPair(76+charge, sqsFeatures[charge][1] * oneOverTotalIntensity); // SQS_INTEN_IN_PAIRS_CHARGE...

		for (int charge=1; charge<=maxChargeForLoop; charge++)
			sample.addPair(80+charge, sqsFeatures[charge][3]); // SQS_NUM_C2_PAIRS_CHARGE...

		for (int charge=1; charge<=maxChargeForLoop; charge++)
			sample.addPair(84+charge, sqsFeatures[charge][4] * oneOverTotalIntensity); // SQS_INTEN_IN_C2_PAIRS_CHARGE...

		for (int charge=1; charge<=maxChargeForLoop; charge++)
			sample.addPair(88+charge, sqsFeatures[charge][0]-sqsFeatures[charge][3]); // SQS_DIFF_NUM_PAIRS_CHARGE...

		for (int charge=1; charge<=maxChargeForLoop; charge++)
			sample.addPair(92+charge, (sqsFeatures[charge][1]+sqsFeatures[charge][4]) *
										oneOverTotalIntensity); // SQS_SUM_INTENSITY_IN_PAIRS_CHARGE1
	}

	// Relative intensity features
	if (1)
	{
		intensity_t maxInten=0.0;
		for (size_t i=0; i<numPeaks; i++)
			if (peaks[i].intensity>maxInten)
				maxInten = peaks[i].intensity;
		
		const intensity_t thresh025=0.25*maxInten, thresh010=0.1*maxInten, thresh001=0.01*maxInten;
		const float total001=totalIntensity*0.01, total005=totalIntensity*0.05;
		int num025=0, num010=0, num001=0;
		int numTotal001=0, numTotal005=0;
		for (size_t i=0; i<numPeaks; i++)
		{
			const intensity_t inten = peaks[i].intensity;
			if (inten>total001)
			{
				numTotal001++;
				if (inten>total005)
					numTotal005++;
			}
			if (inten<thresh001)
				continue;
			num001++;
			if (inten<thresh010)
				continue;
			num010++;
			if (inten<thresh025)
				continue;
			num025++;
		}

		sample.addPair(97, currentSpectrumNumStrongPeaks_ / currentSpectrumTotalIntensity_); 
		sample.addPair(98, num025); // SQS_PROP_PEAKS_025_REL
		sample.addPair(99, num010); // SQS_PROP_PEAKS_010_REL
		sample.addPair(100, num001); // SQS_PROP_PEAKS_001_REL
		sample.addPair(101, numTotal001); // SQS_PROP_PEAKS_1%_TOTAL
		sample.addPair(102, numTotal005); // SQS_PROP_PEAKS_5%_TOTAL
	}

		// detla mass features
	if (1 && numPeaks>6)
	{
		static vector< IntenMassPair > sortedPeaks;
		static vector< mass_t > deltas;
		sortedPeaks.resize(numPeaks);
		sortedPeaks[0].inten = peaks[0].intensity;
		sortedPeaks[0].mass  = peaks[0].mass;

		deltas.clear();
		deltas.reserve(numPeaks-1);
		MeanSdStats msd;
		for (size_t i=1; i<numPeaks; i++)
		{
			const mass_t delta = peaks[i].mass - peaks[i-1].mass;
			deltas.push_back(delta);
			msd.addWX(1.0, delta);
			sortedPeaks[i].inten = peaks[i].intensity;
			sortedPeaks[i].mass  = peaks[i].mass;
		}
		sort(sortedPeaks.begin(), sortedPeaks.end());

		intensity_t sumTop = 0.0;
		const size_t numTop = numPeaks/3;
		for (size_t i=0; i<numTop; i++)
			sumTop += sortedPeaks[i].inten;

		sample.addPair(103, sumTop/totalIntensity); // SQS_PROP_INTEN_THIRD_TOP
		
		double mean=0, sd=0;
		msd.calcMeanAndSd(mean, sd);
		sample.addPair(104, mean); // SQS_AVG_PEAK_DELTA_MASS
		sample.addPair(105, sd);  // SQS_STD_PEAK_DELTA_MASS

		sort(deltas.begin(),deltas.end());
		sample.addPair(106, deltas[deltas.size()/2]); // SQS_MED_PEAK_DELTA_MASS
		// assume 1/6 peaks are strong
		if (numPeaks>=25)
		{
			static vector<mass_t> sortedStrongMasses;
			const size_t numStrong = numPeaks/6;
			sortedStrongMasses.resize(numStrong);
			for (size_t i=0; i<numStrong; i++)
				sortedStrongMasses[i]=sortedPeaks[i].mass; // take the mass

			sort(sortedStrongMasses.begin(), sortedStrongMasses.end());
			MeanSdStats msd;
			for (size_t i=1; i<numStrong; i++)
				msd.addWX(1.0, sortedStrongMasses[i] - sortedStrongMasses[i-1]);

			double mean=0, sd=0;
			msd.calcMeanAndSd(mean, sd);
			sample.addPair(107, mean); // SQS_AVG_STRONG_PEAK_DELTA_MASS
			sample.addPair(108, sd);  // SQS_STD_STRONG_PEAK_DELTA_MASS
		}
	}

	assert( sample.checkConsistency());
}








void calculateBackgroundStatistics(const mass_t single_charge_pair_sum, // the sum of b+y or c+z
						  const Config *config,
						  const PeakList& pl,
						  const vector<bool>& strong_inds,
						  const vector<float>& iso_levels,
						  const vector<bool>& strict_iso_inds,
						  PMCRankStats& pmc_stats_total)
{
	const mass_t tolerance = config->getTolerance();
	const mass_t offsets[]={-22.0,-10.0,-8.5,8.5,12.0,22.5};
	const int numOffsets = sizeof(offsets)/sizeof(mass_t);
	const float one_over_num_offsets = 1.0 / (float)numOffsets;

	int i;
	for (i=0; i<numOffsets; i++)
	{
		static PMCRankStats pmcStats; 
		
		calculatePmcRankStatisticsForMass(pl, single_charge_pair_sum+offsets[i], 
			1.5*tolerance, iso_levels, strict_iso_inds, strong_inds, pmcStats);

		if (i==0)
		{
			pmc_stats_total = pmcStats;
		}
		else
		{
			pmc_stats_total.numCharge2FragmentPairs        += pmcStats.numCharge2FragmentPairs;
			pmc_stats_total.numFragmentPairs           += pmcStats.numFragmentPairs;
			pmc_stats_total.numStrongCharge2FragmentPairs += pmcStats.numStrongCharge2FragmentPairs;
			pmc_stats_total.numStrongFragmentPairs	 += pmcStats.numStrongFragmentPairs;

			pmc_stats_total.intensityInFragmentPairs		 += pmcStats.intensityInFragmentPairs;
			pmc_stats_total.intensityInStrongPairs		 += pmcStats.intensityInStrongPairs;
			pmc_stats_total.intensityInCharge2Pairs			 += pmcStats.intensityInCharge2Pairs;
			pmc_stats_total.intensityInStrongCharge2Pairs	 += pmcStats.intensityInStrongCharge2Pairs;

			
		}
	}

	pmc_stats_total.numCharge2FragmentPairs		 *= one_over_num_offsets;
	pmc_stats_total.numFragmentPairs			 *= one_over_num_offsets;
	pmc_stats_total.numStrongCharge2FragmentPairs *= one_over_num_offsets;
	pmc_stats_total.numStrongFragmentPairs	 *= one_over_num_offsets;

	pmc_stats_total.intensityInFragmentPairs		 *= one_over_num_offsets;
	pmc_stats_total.intensityInStrongPairs		 *= one_over_num_offsets;
	pmc_stats_total.intensityInCharge2Pairs			 *= one_over_num_offsets;
	pmc_stats_total.intensityInStrongCharge2Pairs	 *= one_over_num_offsets;

}





/************************************************************************
Looks at the best m/z value for each charge and returns
the number of fragment pairs, strong fragment pairs, 
charge 2 fragment pairs, and strong charge 2 fragment pairs.
The selection is done according to what m/z has the maximal sum of pairs+
strong pairs and also what m/z has the maximal sum of c2 pairs + strong c2
pairs
*************************************************************************/
void PMCSQS_Scorer::getSqsFeaturesFromPmcTabels(const PeakList& pl, 
									            vector< vector<float> >& sqsFeatures) const
{
	const mass_t originalMOverZ = pl.getHeader()->getMOverZ();
	
	if (sqsFeatures.size() != maximalChargeWithModels_ +1)
		sqsFeatures.resize(maximalChargeWithModels_+1);
	sqsFeatures[0].clear();
	sqsFeatures[0].resize(6,0.0);

	for (int charge=1; charge<=maximalChargeWithModels_; charge++)
	{
		// find entry which has the maximal number of strong pairs, while maining the minimial
		// m/z distance shift from the original
		// for selecting the charge1 pair statistics, look at the best objective
		// with the charge 1 pairs. Similarly look at charge2 pairs for c2 objective
		
		float maxObjectiveValue1=0.0, maxObjectiveValue2=0.0;
		int   bestTableIndex1=0, bestTableIndex2=0;
		mass_t mzDiff1 = MAX_FLOAT, mzDiff2 = MAX_FLOAT;

		for (size_t i=1; i<currentSpectrumPmcTables_[charge].size(); i++)
		{
			const PMCRankStats& stats = currentSpectrumPmcTables_[charge][i];

			const float currentObjectiveValue1 = stats.numStrongFragmentPairs +
					stats.numStrongFragmentPairs + stats.numFragmentPairs;

			const float currentObjectiveValue2 = stats.numCharge2FragmentPairs +
					stats.numStrongCharge2FragmentPairs + stats.numStrongCharge2FragmentPairs;
			
			if (currentObjectiveValue1 >= maxObjectiveValue1)
			{
				const mass_t distance = fabs(currentSpectrumPmcTables_[charge][i].m_over_z - originalMOverZ);
				if (currentObjectiveValue1 > maxObjectiveValue1 || distance>= mzDiff1)
				{
					maxObjectiveValue1 = currentObjectiveValue1;
					mzDiff1			  = distance;
					bestTableIndex1    = i;
				}
			}

			if (currentObjectiveValue2 >= maxObjectiveValue2)
			{
				const mass_t distance = fabs(currentSpectrumPmcTables_[charge][i].m_over_z - originalMOverZ);
				if (currentObjectiveValue2 > maxObjectiveValue2 || distance>= mzDiff1)
				{
					maxObjectiveValue2 = currentObjectiveValue2;
					mzDiff2			   = distance;
					bestTableIndex2    = i;
				}
			}
		}

		sqsFeatures[charge].resize(6);
		sqsFeatures[charge][0]=currentSpectrumPmcTables_[charge][bestTableIndex1].numFragmentPairs;
		sqsFeatures[charge][1]=currentSpectrumPmcTables_[charge][bestTableIndex1].intensityInFragmentPairs;
		sqsFeatures[charge][2]=currentSpectrumPmcTables_[charge][bestTableIndex1].intensityInStrongPairs;

		sqsFeatures[charge][3]=currentSpectrumPmcTables_[charge][bestTableIndex2].numCharge2FragmentPairs;
		sqsFeatures[charge][4]=currentSpectrumPmcTables_[charge][bestTableIndex2].intensityInCharge2Pairs;
		sqsFeatures[charge][5]=currentSpectrumPmcTables_[charge][bestTableIndex2].intensityInStrongCharge2Pairs;
		

	/*	cout << endl << charge << "\t" << currentSpectrumPmcTables_[charge].size() << "\t" <<
			sqsFeatures[charge][0] << "\t" << sqsFeatures[charge][1] << "\t" <<
			sqsFeatures[charge][2] << "\t" << sqsFeatures[charge][3] << endl;*/

	}
}






/*****************************************************************
Creates for each input file an mgf file that holds the spectra
that passed quality filtering does not correct PM and charge in the 
mgf files.
******************************************************************
void PMCSQS_Scorer::output_filtered_spectra_to_mgfs(
									 Config *config,
									 const vector<string>& files,
									 char *out_dir,
									 float filter_prob, 
									 int& total_num_written, 
									 int& total_num_read)
{
	total_num_read = 0;
	total_num_read = 0;
	int f;
	for (f=0; f<files.size(); f++)
	{
		const char *spectra_file = files[f].c_str();
		FileManager fm;
		FileSet fs;
		BasicSpecReader bsr;

		string fname, mgf_name, map_name;
		getFileNameWithoutExtension(files[f].c_str(), fname);

		mgf_name = string(out_dir) + "/" + fname + "_fil.mgf";
		map_name = string(out_dir) + "/" + fname + "_map.txt";

		///////////////////////////////////////////////
		// Quick read, get all pointers to begining of spectra
		if (get_file_extension_type(files[f]) != MZXML)
		{
			fm.init_from_file(config,spectra_file);
		}
		else
			fm.init_and_read_single_mzXML(config,spectra_file,f);

		fs.select_all_files(fm);

		const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
		int sc;
		int  num_spec_written=0;
		bool first=true;
		ofstream out_stream, map_stream;
		for (sc=0; sc<all_ssf.size(); sc++)
		{
			static vector<QCPeak> peaks;
			SingleSpectrumFile *ssf = all_ssf[sc];
			
			if (peaks.size()<ssf->num_peaks)
			{
				int new_size = ssf->num_peaks*2;
				if (new_size<2500)
					new_size=2500;
				peaks.resize(new_size);
			}

			// read without processing peaks
			const int num_peaks = bsr.read_basic_spec(config,fm,ssf,&peaks[0],false,true);
			ssf->file_idx = f;

			BasicSpectrum bs;
			bs.peaks = &peaks[0];
			bs.num_peaks = num_peaks;
			bs.ssf = all_ssf[sc];

			int max_charge;
			float prob=0.0; // = get_sqs_for_spectrum(config,bs,&max_charge);
			if (prob<filter_prob)
			{
				continue;
			}

			if (first)
			{
				out_stream.open(mgf_name.c_str(),ios::out);
				map_stream.open(map_name.c_str(),ios::out);

				cout << "Filtering spectra to minumum quality score: " << filter_prob << endl;
				cout << "Writing spectra info to:" << endl;
				cout << mgf_name << endl << map_name << endl;

				if (! out_stream.is_open() || ! out_stream.good())
				{
					cout << "Error: couldn\'t open for out mgf stream for writing: " <<
						endl << mgf_name << endl;
					exit(1);
				}
				first = false;
			}

			char single_name[64];
			sprintf(single_name,"%d:%d",f,bs.ssf->get_scan());
			bs.ssf->single_name = single_name;
			bs.output_to_mgf(out_stream,config);
			if (prob>1.0)
				prob=1.0;
			map_stream << num_spec_written++ << "\t" << all_ssf[sc]->get_scan() << "\t" << fixed << prob << endl;

			
		}
		out_stream.close();
		map_stream.close();

		total_num_read+= all_ssf.size();
		total_num_written += num_spec_written;
	}
}*/


	
