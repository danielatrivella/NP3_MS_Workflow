#include "PMCSQS.h"
#include "SpectraList.h"
#include "PepNovo_auxfun.h"
#include "../MLlib/MlRankBoost.h"


bool PMCSQS_Scorer::read_pmc_rank_models(const Config *_config, const char *file_name)
{
	config_ = _config;

	string path = config_->get_resource_dir() + "/" + file_name;
	ifstream in_stream(path.c_str(),ios::in);
	if (! in_stream.good())
	{
		cout << "Warning: couldn't open pmc rank model for reading: " << path << endl;
		return false;
	}


	char buff[512];
	int numCharges=-1;

	in_stream.getline(buff,256);
	istringstream iss1(buff);

	fragmentPairSumOffset=NEG_INF;
	bin_increment=NEG_INF;
	iss1 >> bin_increment >> this->fragmentPairSumOffset;
	if (fragmentPairSumOffset==NEG_INF || bin_increment == NEG_INF)
	{
		cout << "Error in pmc model file!" << endl;
		exit(1);
	}

	in_stream.getline(buff,256);
	istringstream iss(buff);

	iss >> numCharges;
	maximalChargeWithModels_=numCharges-1;
	
	pmc_rank_models.resize(numCharges);
	pmcMassThresholds_.resize(numCharges);
	pmcMzBiases_.resize(numCharges);


	int i;
	for (i=0; i<numCharges; i++)
	{
		in_stream.getline(buff,256);
		istringstream iss(buff);
		int num_threshes=0;
		iss >> num_threshes;
		
		pmcMassThresholds_[i].resize(num_threshes,NEG_INF);
		int j;
		for (j=0; j<num_threshes; j++)
			iss >> pmcMassThresholds_[i][j];
	}

	for (i=0; i<numCharges; i++)
	{
		in_stream.getline(buff,256);
		istringstream iss(buff);
		int num_biases=0;
		iss >> num_biases;
		
		pmcMzBiases_[i].resize(num_biases,NEG_INF);
		int j;
		for (j=0; j<num_biases; j++)
			iss >> pmcMzBiases_[i][j];
	}
	
	// read Boost models
	for (i=0; i<numCharges; i++)
	{
		in_stream.getline(buff,256);
		istringstream iss(buff);

		int num_models=-1;
		iss >> num_models;

		if (num_models<0)
		{
			cout << "Error: bad parsing of PMC model file!" << endl;
			exit(0);
		}
		pmc_rank_models[i].resize(num_models,NULL);

		int j;
		for (j=0; j<num_models; j++)
		{
			pmc_rank_models[i][j]=new RankBoostModel;
			pmc_rank_models[i][j]->read_rankboost_model(in_stream);
		}
		
	}
	in_stream.close();


	indInitializedPmc_ = true;
	return true;
}


void PMCSQS_Scorer::write_pmc_rank_models(const char *path) const
{
	ofstream out_stream(path,ios::out);
	if (! out_stream.good())
	{
		cout << "Error: couldn't open pmc model for writing: " << path << endl;
		exit(1);
	}

	out_stream << this->bin_increment << " " << this->fragmentPairSumOffset << endl;
	out_stream << this->pmcMassThresholds_.size() << endl;
	out_stream << setprecision(3) << fixed;

	int i;
	for (i=0; i<this->pmcMassThresholds_.size(); i++)
	{
		out_stream << pmcMassThresholds_[i].size();
		int j;
		for (j=0; j<pmcMassThresholds_[i].size(); j++)
			out_stream << " " << pmcMassThresholds_[i][j];
		out_stream << endl;
	}

	
	for (i=0; i<this->pmcMzBiases_.size(); i++)
	{
		out_stream << pmcMzBiases_[i].size();
		int j;
		for (j=0; j<pmcMzBiases_[i].size(); j++)
			out_stream << " " << pmcMzBiases_[i][j];
		out_stream << endl;
	}

	// write RankBoost models
	for (i=0; i<pmc_rank_models.size(); i++)
	{
		int j;
		
		if (pmc_rank_models[i].size()==1 && ! pmc_rank_models[i][0])
		{
			out_stream << 0 << endl;
			continue;
		}

		out_stream << pmc_rank_models[i].size() << endl;
		for (j=0; j<pmc_rank_models[i].size(); j++)
		{
			if (pmc_rank_models[i][j])
			{
				pmc_rank_models[i][j]->write_rankboost_model(out_stream,true);
			}
			else
			{
				cout << "Error: non intialized rank pmc model!" << endl;
				exit(1);
			}
		}
	}
	
	out_stream.close();
}


void PMCSQS_Scorer::setPmcrMassThresholds(int option)
{
	if (option==0)
	{
		pmcMassThresholds_ = config_->get_size_thresholds();
		for (size_t i=0; i<pmcMassThresholds_.size(); i++)
		{
			if (pmcMassThresholds_[i].size() > 0 && 
				pmcMassThresholds_[i].back() > 9999.0)
				pmcMassThresholds_[i].pop_back();
		}
	}

	if (option==1)
	{
		pmcMassThresholds_.resize(4);
		pmcMassThresholds_[1].push_back(1150.0);
		pmcMassThresholds_[1].push_back(1400.0);
 
		pmcMassThresholds_[2].push_back(1100.0);
		pmcMassThresholds_[2].push_back(1300.0);
		pmcMassThresholds_[2].push_back(1600.0);
		pmcMassThresholds_[2].push_back(1900.0);
		pmcMassThresholds_[2].push_back(2400.0);

		pmcMassThresholds_[3].push_back(1950.0);
		pmcMassThresholds_[3].push_back(2450.0);
		pmcMassThresholds_[3].push_back(3000.0);
	}
}





void PMCSQS_Scorer::trainPmcRankModels(const Config* config,
									   const SpectraAggregator& sa, 
									   int specificCharge,
									   bool overwrite)
{
	const vector<int>& spectraCounts = sa.getSpectraCountsPerCharge();
	const int maxSpectraPerFile = 40000;
	
	string featureFile = config->get_resource_dir() + "/" + config->getModelName() + "_PMC_fet.txt";
	MlFeatureSet fs;
	fs.readFeatureFile(featureFile.c_str());

	// try and read existing pmc model, otherwise init a new one
	string pmcPath = config->get_resource_dir() + "/" + config->getModelName() + "_PMC.txt";
	ifstream modelStream(pmcPath.c_str());	

	if (! modelStream.is_open() && modelStream.good())
	{

		modelStream.close();
		string pmcrName = config->getModelName() + "_PMC.txt";
	//	readPmcrModels(config, pmcrName.c_str());
	}
	else
	{
		setPmcrMassThresholds();
		setFragmentPairSumOffset(MASS_PROTON); // b+y - PM+19
		setBinMassIncrement(0.1);
		
		pmcRankBoostModels_.resize(pmcMassThresholds_.size());
		pmcMzBiases_.resize(pmcMassThresholds_.size());
	}
	
	maximalChargeWithModels_= pmcMassThresholds_.size()-1;
	const double proportionForTraining = 0.5;


	// It is assumed that the mass thresholds were set according to the training data
	// (this is done manually with values encoded in the set_mass_threhsolds function)
	for (size_t charge=1; charge<=maximalChargeWithModels_; charge++)
	{
		if (specificCharge>0 && charge != specificCharge)
			continue;

		const size_t numSizes = pmcMassThresholds_[charge].size();
		pmcRankBoostModels_[charge].resize(numSizes+1, NULL);
		pmcMzBiases_[charge].resize(numSizes+1, 0.0);

		for (size_t sizeIndex=0; sizeIndex<=numSizes; sizeIndex++)
		{
			mass_t minMz = 0.0;
			mass_t maxMz = MAX_FLOAT;

			if (sizeIndex>0)
				minMz = pmcMassThresholds_[charge][sizeIndex-1];

			if (sizeIndex<numSizes)
				maxMz = pmcMassThresholds_[charge][sizeIndex];

			minMz /= static_cast<mass_t>(charge);
			maxMz /= static_cast<mass_t>(charge);

			cout << "\nTraining PMC rank model for charge " << charge << " size " << sizeIndex << endl;
			cout << "min m/z: " << minMz << endl;
			cout << "max m/z: " << maxMz << endl;
			if (pmcRankBoostModels_[charge][sizeIndex] && ! overwrite)
			{
				cout << endl << "Already trained..." << endl;
				continue;
			}

			vector<const SingleSpectrumHeader*> testHeaders;
			MlRankBoostDataset trainDS, testDS; 

			// check for existing training/test data files use them if they exist
			ostringstream oss ;
			oss << config->getTrainingDir() << "/" << config->getModelName() << "_PMC_" << charge << "_" << sizeIndex;
			string trainFile = oss.str() + "_tr.txt";
			string testFile  = oss.str() + "_ts.txt";
			
			if (checkIfFileExists(trainFile.c_str()) && checkIfFileExists(testFile.c_str()))
			{
				const size_t numTrainRead = trainDS.readRankBoostDataSet(trainFile.c_str());
				const size_t numTestRead  = testDS.readRankBoostDataSet(testFile.c_str());

				assert(trainDS.getPhi().size()>0);
				assert(testDS.getPhi().size()>0);
				cout << "Using existing data files:" << endl;
				cout << "Read " << numTrainRead << " training samples from " << trainFile << " (including "
					 << trainDS.getPhi().size() << " pairs for the feedback function)" << endl;
				cout << "Read " << numTestRead << " test samples from " << testFile << " (including "
					 << testDS.getPhi().size() << " pairs for the feedback function)" << endl << endl;
			}
			else
			{
				// these ranges are given according to pm_with_19
				// so files should be selected through select_files and not
				// select_file_in_mz_range
				SpectraList sl(sa);
				sl.selectHeaders(minMz, maxMz, charge, charge);
			

				if (sl.getNumHeaders()< MIN_SPECTRA_FOR_PMCSQS_MODEL)
				{
					cout << "Not enough spectra to train this model, found only " << sl.getNumHeaders() << ", skipping..." << endl;
					continue;
				}
				sl.randomlyReduceListToSize(maxSpectraPerFile);

				cout << "PM m/z range " << minMz << "-" << maxMz << "  using " << sl.getNumHeaders() << " spectra." << endl;

				size_t numTrainingGroups=0;
				size_t numTestingGroups=0;
				
				const size_t numSamples = sl.getNumHeaders();
				
				// first find the bias in number of bins between the true m/z bin and
				// the optimal m/z bin
				vector<bool> skippedIndexes(numSamples, false);
				size_t skippedBadMz=0;
				mass_t totalBias=0;
				for (size_t i=0; i<numSamples; i++)
				{
					const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
					if (header->getPeptideStr().length()<=0)
						continue;

					PeakList pl;
					pl.readPeaksToLocalAllocation(sa, header);
					Peptide peptide;
					peptide.parseFromString(config, header->getPeptideStr());
					
					const mass_t trueMz = (peptide.get_mass()+MASS_H2O+ charge*MASS_PROTON)/static_cast<mass_t>(charge);

					if (fabs(trueMz - header->getMOverZ()) > 2.5)
					{
						skippedBadMz++;
						skippedIndexes[i]=true;
						continue;
					} 

					initializeForCurrentSpectrum(config, pl);
					calculateCurrentSpectrumPmcValues(pl, bin_increment);
		
					// find the trueMzBinIdx
					const vector<PMCRankStats>& pmcStats = currentSpectrumPmcTables_[charge];
					size_t trueMzBinIdx=0;
					while (trueMzBinIdx<pmcStats.size() && pmcStats[trueMzBinIdx].m_over_z<trueMz)
						trueMzBinIdx++;

					if (trueMzBinIdx == pmcStats.size())
						trueMzBinIdx--;

					if (trueMzBinIdx>0 && pmcStats[trueMzBinIdx].m_over_z-trueMz>trueMz - pmcStats[trueMzBinIdx-1].m_over_z)
						trueMzBinIdx--;

					const size_t optBinIdx = findOptimalBinIdx(trueMzBinIdx, charge);
					if (optBinIdx ==0 || optBinIdx >= pmcStats.size()-1)
					{
						skippedBadMz++;
						skippedIndexes[i]=true;
						continue;
					}

					totalBias += (pmcStats[optBinIdx].m_over_z - pmcStats[trueMzBinIdx].m_over_z);

					if (fabs(pmcStats[optBinIdx].m_over_z - pmcStats[trueMzBinIdx].m_over_z)>4.0)
					{
						cout << "opt bin: " << optBinIdx << " (" << pmcStats[optBinIdx].m_over_z << ")  ";
						cout << "tru bin: " << trueMzBinIdx << " ("<< pmcStats[trueMzBinIdx].m_over_z << ")" << endl;
					}
				} 
		
				assert(numSamples>skippedBadMz);
				mass_t mzBias = totalBias / static_cast<mass_t>(numSamples-skippedBadMz);
				pmcMzBiases_[charge][sizeIndex]=mzBias;

				cout << "m/z bias: " << setprecision(4) << mzBias << endl;
				cout << "skipped " << skippedBadMz << "/" << numSamples <<
					"  because of m/z more than 2.5 away from observed..." << endl; 

				vector<MlSample> spectrumRankSamples;
				for (size_t i=0; i<numSamples; i++)
				{
					if (skippedIndexes[i])
						continue;

					const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
					PeakList pl;

					pl.readPeaksToLocalAllocation(sa, header);
					if (header->getPeptideStr().length() == 0)
						continue;

					Peptide peptide;
					peptide.parseFromString(config, header->getPeptideStr());
					
					const mass_t trueMz = (peptide.get_mass()+MASS_H2O+ charge*MASS_PROTON)/static_cast<mass_t>(charge);

					initializeForCurrentSpectrum(config, pl);
					calculateCurrentSpectrumPmcValues(pl, bin_increment);

					// find the trueMzBinIdx
					const vector<PMCRankStats>& pmcStats = currentSpectrumPmcTables_[charge];
					int trueMzBinIdx=0;
					while (trueMzBinIdx<pmcStats.size() && pmcStats[trueMzBinIdx].m_over_z<trueMz)
						trueMzBinIdx++;

					if (trueMzBinIdx == pmcStats.size())
						trueMzBinIdx--;

					if (trueMzBinIdx>0 && pmcStats[trueMzBinIdx].m_over_z-trueMz>trueMz-pmcStats[trueMzBinIdx-1].m_over_z)
						trueMzBinIdx--;

					const size_t optBinIdx = findOptimalBinIdx(trueMzBinIdx, charge);

					fillRankboostPmcSamples(pl, charge, spectrumRankSamples);

					// select samples and add them to pmc_ds
					size_t goodIndex=MAX_SIZE_T;
					vector<size_t> badIndexes;
					selectTrainingSampleIndexes(charge, spectrumRankSamples, pl, goodIndex, badIndexes);
					assert(goodIndex<MAX_SIZE_T);

					const bool indAddToTraining = (myRandom()<proportionForTraining);
					size_t groupIndex;
					if (indAddToTraining)
					{
						groupIndex= numTrainingGroups++;	
					}
					else
					{
						groupIndex= numTestingGroups++;
						testHeaders.push_back(header);
					}
					
					
					MlRankBoostDataset& ds = (indAddToTraining ? trainDS : testDS);

					const int positiveIndex  = ds.getNumSamples();
					spectrumRankSamples[goodIndex].groupIndex  = groupIndex;
					spectrumRankSamples[goodIndex].rankInGroup = 0;
					spectrumRankSamples[goodIndex].label = 0;
					ds.addSample(spectrumRankSamples[goodIndex]);
		
					for (size_t j=0; j<badIndexes.size(); j++)
					{
						const size_t badIndex = badIndexes[j];
						assert( badIndex < spectrumRankSamples.size() );
			
						spectrumRankSamples[badIndex].groupIndex  = groupIndex;
						spectrumRankSamples[badIndex].rankInGroup = 1;
						spectrumRankSamples[badIndex].label = 1;
						ds.addPairToPhiVector(ds.getNumSamples(),positiveIndex);
						ds.addSample(spectrumRankSamples[badIndex]);
					}				   
				}
				
				trainDS.writeRankBoostDataSet(trainFile.c_str());
				testDS.writeRankBoostDataSet(testFile.c_str());
			}
				
			if (pmcRankBoostModels_[charge][sizeIndex])
				delete pmcRankBoostModels_[charge][sizeIndex];

			pmcRankBoostModels_[charge][sizeIndex] = new MlModel;

			MlTrainingContainer mtc; // all training options go here

			mtc.setTrainingAlgorithm(SMT_RANKBOOST_ADA);
			mtc.setTrainingDataSet(&trainDS);
			mtc.setTestDataSet(&testDS);
			mtc.setMaxNumIterations(10000);
			mtc.setOutputDir("Training");
			mtc.setOutputName(oss.str().c_str());

			pmcRankBoostModels_[charge][sizeIndex]->createModel(mtc);
		}
	}

//	string path;
//	path = config->get_resource_dir() + "/" + config->getModelName() + "_PMC.txt";
//	this->write_pmc_rank_models(path.c_str());
//	indInitializedPmc_ = true; */
}







struct offset_pair {
	offset_pair() : offset(POS_INF), inten_sum(0) {};
	offset_pair(mass_t off,float inten) : offset(off), inten_sum(inten) {};
	mass_t offset;
	float inten_sum;
};

bool cmp_offset_pair_offset (const offset_pair& a, const offset_pair& b)
{
	return (a.offset<b.offset);
}

bool cmp_offset_pair_inten (const offset_pair& a, const offset_pair& b)
{
	return (a.inten_sum>b.inten_sum);
}


float calc_mean_abs_offset(const vector<float>& offsets_by_inten)
{
	const float missing_pair_offset = 0.5;
	const int   num_offsets         = 3;

	if (offsets_by_inten.size()==0)
		return 1000;

	float abs_off=0;
	int i;
	for (i=0; i<num_offsets && i<offsets_by_inten.size(); i++)
		abs_off+=fabs(offsets_by_inten[i]);

	abs_off += (3-i)*missing_pair_offset;
	
	return (abs_off/num_offsets);
}


void calculatePmcRankStatisticsForMass(const PeakList& pl, 
									   mass_t SingleChargeFragmentPairSum,
									   mass_t tolerance, 
									   const vector<float>& isotopicLevels,
									   const vector<bool>& strictIsotopicIndicators,
									   const vector<bool>& strongPeakInds,
									   PMCRankStats& statistics)
{
	const mass_t minSingleSum = SingleChargeFragmentPairSum - tolerance;
	const mass_t maxSingleSum = SingleChargeFragmentPairSum + tolerance;

	const mass_t minDoubleSum = minSingleSum + MASS_PROTON;
	const mass_t maxDoubleSum = maxSingleSum + MASS_PROTON;
	const mass_t doubleChargeFragmentPairSum = SingleChargeFragmentPairSum + MASS_PROTON;

	const mass_t minSingleSumWithH2OLoss = minSingleSum - MASS_H2O;
	const mass_t maxSingleSumWithH2OLoss = maxSingleSum - MASS_H2O;
	const mass_t single_charge_pair_h2o_sum = SingleChargeFragmentPairSum - MASS_H2O;

	const mass_t min_double_h2o_sum = minDoubleSum - MASS_H2O;
	const mass_t max_double_h2o_sum = maxDoubleSum - MASS_H2O;
	const mass_t double_charge_pair_h2o_sum = doubleChargeFragmentPairSum - MASS_H2O;

	static vector<offset_pair> by_pairs,  strong_pairs;
	static vector<offset_pair> c2_pairs,  strong_c2_pairs;
	static vector<offset_pair> h2o_pairs, c2_h2o_pairs;

	by_pairs.clear();
	strong_pairs.clear();
	c2_pairs.clear();
	strong_c2_pairs.clear();
	h2o_pairs.clear();
	c2_h2o_pairs.clear();

	statistics.clear();

	const Peak* const peaks = pl.getPeaks();
	const int numPeaks = pl.getNumPeaks();
	int forwardIndex = -1;
	int backIndex = numPeaks-1;


	// find pairs of b/y
	while (forwardIndex<backIndex)
	{
		if (isotopicLevels[++forwardIndex]>0)
			continue;

		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>maxSingleSum)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > minSingleSum)
		{
			if (isotopicLevels[backIndex]>0)
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - SingleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			by_pairs.push_back(offset_pair(offset,inten_sum));
			statistics.intensityInFragmentPairs += inten_sum;

			if (strongPeakInds[forwardIndex] || strongPeakInds[backIndex])
			{
				strong_pairs.push_back(offset_pair(offset,inten_sum));
				statistics.intensityInStrongPairs += inten_sum;
			}
		}
	}

	// find pairs b/y2
	forwardIndex = -1;
	backIndex = numPeaks-1;

	const int last_idx = backIndex;
	while (forwardIndex<last_idx && backIndex>=0)
	{
		if (isotopicLevels[++forwardIndex]>0)
			continue;
			
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (sum>maxDoubleSum)
		{
			if (--backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > minDoubleSum)
		{
			if (isotopicLevels[backIndex]>0)
				continue;

			const mass_t offset = fabs(sum - doubleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_pairs.push_back(offset_pair(offset,inten_sum));
			statistics.intensityInCharge2Pairs += inten_sum;

			if (strongPeakInds[forwardIndex] || strongPeakInds[backIndex])
			{
				strong_c2_pairs.push_back(offset_pair(offset,inten_sum));
				statistics.intensityInStrongCharge2Pairs = inten_sum;
			}
		}
	}

	// find pairs of b/y-H2O
	forwardIndex = -1;
	backIndex = numPeaks-1;

	while (forwardIndex<backIndex)
	{
		if (isotopicLevels[++forwardIndex]>0)
			continue;

		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>maxSingleSumWithH2OLoss)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > minSingleSumWithH2OLoss)
		{
			if (isotopicLevels[backIndex]>0)
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - single_charge_pair_h2o_sum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			h2o_pairs.push_back(offset_pair(offset,inten_sum));
			statistics.intensityInH2oLossPairs += inten_sum;
		}
	}

	// find pairs b/y2 - H2O
	forwardIndex = -1;
	backIndex = numPeaks-1;

	while (forwardIndex<last_idx && backIndex>=0)
	{
		if (isotopicLevels[++forwardIndex]>0)
			continue;
			
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (sum>max_double_h2o_sum)
		{
			if (--backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > min_double_h2o_sum)
		{
			if (isotopicLevels[backIndex]>0)
				continue;

			const mass_t offset = fabs(sum - double_charge_pair_h2o_sum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_h2o_pairs.push_back(offset_pair(offset,inten_sum));
			statistics.intensityInCharge2H2oLossPairs += inten_sum;
		}
	}

	statistics.numFragmentPairs = by_pairs.size();
	statistics.numStrongFragmentPairs = strong_pairs.size();
	statistics.numCharge2FragmentPairs = c2_pairs.size();
	statistics.numStrongCharge2FragmentPairs = strong_c2_pairs.size();
	statistics.numH2oLossFragmentPairs = h2o_pairs.size();
	statistics.numCharge2H2oLossFragmentPairs = c2_h2o_pairs.size();

	vector<float>& offsetPairsOderedByIntensity = statistics.offsetPairsOderedByIntensity;
	sort(by_pairs.begin(),by_pairs.end(),cmp_offset_pair_inten);
	offsetPairsOderedByIntensity.resize(by_pairs.size());
	for (size_t i=0; i<by_pairs.size(); i++)
		offsetPairsOderedByIntensity[i]=by_pairs[i].offset;
	statistics.meanOffsetPairs=calc_mean_abs_offset(offsetPairsOderedByIntensity);

	vector<float>& strong_offsetPairsOderedByIntensity = statistics.strong_offsetPairsOderedByIntensity;
	sort(strong_pairs.begin(),strong_pairs.end(),cmp_offset_pair_inten);
	strong_offsetPairsOderedByIntensity.resize(strong_pairs.size());
	for (size_t i=0; i<strong_pairs.size(); i++)
		strong_offsetPairsOderedByIntensity[i]=strong_pairs[i].offset;
	statistics.meanOffsetStrongPairs=calc_mean_abs_offset(strong_offsetPairsOderedByIntensity);

	vector<float>& charge2offsetPairsOderedByIntensity = statistics.charge2offsetPairsOderedByIntensity;
	sort(c2_pairs.begin(),c2_pairs.end(),cmp_offset_pair_inten);
	charge2offsetPairsOderedByIntensity.resize(c2_pairs.size());
	for (size_t i=0; i<c2_pairs.size(); i++)
		charge2offsetPairsOderedByIntensity[i]=c2_pairs[i].offset;
	statistics.meanOffsetCharge2Pairs=calc_mean_abs_offset(charge2offsetPairsOderedByIntensity);



	// fill in additional iso sum features (look at pairs that sum to expected, expected+1 expected+2)
	// find pairs of b/y

	static vector<offset_pair> pairs0,  pairs1, pairs2;
	static vector<offset_pair> c2_pairs0,  c2_pairs1, c2_pairs2;
	
	pairs0.clear();
	forwardIndex = -1;
	backIndex = numPeaks-1;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
		if (strictIsotopicIndicators[forwardIndex])
			continue;

		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>maxSingleSum)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > minSingleSum)
		{
			if (strictIsotopicIndicators[backIndex])
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - SingleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			pairs0.push_back(offset_pair(offset,inten_sum));
		}
	}

	pairs1.clear();
	forwardIndex = -1;
	backIndex = numPeaks-1;
	const mass_t max1 = maxSingleSum+1.0;
	const mass_t min1 = minSingleSum+1.0;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
	
		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>max1)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > min1)
		{
			if (! (strictIsotopicIndicators[backIndex] || strictIsotopicIndicators[forwardIndex]))
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - SingleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			pairs1.push_back(offset_pair(offset,inten_sum));
		}
	}

	pairs2.clear();
	forwardIndex = -1;
	backIndex = numPeaks-1;
	const mass_t max2 = maxSingleSum+2.0;
	const mass_t min2 = minSingleSum+2.0;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
	
		while (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass>max2)
			backIndex--;

		if (backIndex>=0 && peaks[forwardIndex].mass + peaks[backIndex].mass > min2)
		{
			if (! (strictIsotopicIndicators[backIndex] || strictIsotopicIndicators[forwardIndex]))
				continue;

			const mass_t offset = fabs(peaks[forwardIndex].mass + peaks[backIndex].mass - SingleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
					
			pairs2.push_back(offset_pair(offset,inten_sum));
		}
	}


	c2_pairs0.clear();
	forwardIndex = -1;
	backIndex = numPeaks-1;
	while (forwardIndex<backIndex)
	{
		if (strictIsotopicIndicators[++forwardIndex])
			continue;
			
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (backIndex>=0 && sum>maxDoubleSum)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > minDoubleSum)
		{
			if (strictIsotopicIndicators[backIndex])
				continue;

			const mass_t offset = fabs(sum - doubleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_pairs0.push_back(offset_pair(offset,inten_sum));
		}
	}

	c2_pairs1.clear();
	const mass_t maxc21 = maxDoubleSum + 1.0;
	const mass_t minc21 = minDoubleSum + 1.0;
	forwardIndex = -1;
	backIndex = numPeaks-1;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
	
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (backIndex>=0 && sum>maxc21)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > minc21)
		{
			if (! (strictIsotopicIndicators[backIndex] || strictIsotopicIndicators[forwardIndex]) )
				continue;

			const mass_t offset = fabs(sum - doubleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_pairs1.push_back(offset_pair(offset,inten_sum));
		}
	}


	c2_pairs2.clear();
	const mass_t maxc22 = maxDoubleSum + 2.0;
	const mass_t minc22 = minDoubleSum + 2.0;
	forwardIndex = -1;
	backIndex = numPeaks-1;
	while (forwardIndex<backIndex)
	{
		forwardIndex++;
	
		mass_t sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		while (backIndex>=0 && sum>maxc22)
		{
			backIndex--;
			if (backIndex<0)
				break;
			sum = 2*peaks[forwardIndex].mass + peaks[backIndex].mass;
		}

		if (backIndex>=0 && sum > minc22)
		{
			if (! (strictIsotopicIndicators[backIndex] || strictIsotopicIndicators[forwardIndex]) )
				continue;

			const mass_t offset = fabs(sum - doubleChargeFragmentPairSum);
			const float inten_sum = peaks[forwardIndex].intensity + peaks[backIndex].intensity;
			
			c2_pairs2.push_back(offset_pair(offset,inten_sum));
		}
	}

	// use the first 4 peaks
	size_t maxIdx;

	statistics.intensityInPairsIso0=0;
	statistics.numPairsIso0 = pairs0.size();
	sort(pairs0.begin(),pairs0.end(),cmp_offset_pair_inten);
	maxIdx = (pairs0.size()<4 ? pairs0.size() : 4);
	for (size_t i=0; i<maxIdx; i++)
		statistics.intensityInPairsIso0+=pairs0[i].inten_sum;

	statistics.intensityInPairsIso1=0;
	statistics.numPairsIso1 = pairs1.size();
	sort(pairs1.begin(),pairs1.end(),cmp_offset_pair_inten);
	maxIdx = (pairs1.size()<4 ? pairs1.size() : 4);
	for (size_t i=0; i<maxIdx; i++)
		statistics.intensityInPairsIso1+=pairs1[i].inten_sum;
	
	statistics.intensityInPairsIso2=0;
	statistics.numPairsIso2 = pairs2.size();
	sort(pairs2.begin(),pairs2.end(),cmp_offset_pair_inten);
	maxIdx = (pairs2.size()<4 ? pairs2.size() : 4);
	for (size_t i=0; i<maxIdx; i++)
		statistics.intensityInPairsIso2+=pairs2[i].inten_sum;
	
	statistics.intensityInCharge2PairsIso0=0;
	statistics.numCharge2PairsIso0 = c2_pairs0.size();
	sort(c2_pairs0.begin(),c2_pairs0.end(),cmp_offset_pair_inten);
	maxIdx = (c2_pairs0.size()<4 ? c2_pairs0.size() : 4);
	for (size_t i=0; i<maxIdx; i++)
		statistics.intensityInCharge2PairsIso0+=c2_pairs0[i].inten_sum;
	
	statistics.intensityInCharge2PairsIso1=0;
	statistics.numCharge2PairsIso1 = c2_pairs1.size();
	sort(c2_pairs1.begin(),c2_pairs1.end(),cmp_offset_pair_inten);
	maxIdx = (c2_pairs1.size()<4 ? c2_pairs1.size() : 4);
	for (size_t i=0; i<maxIdx; i++)
		statistics.intensityInCharge2PairsIso1+=c2_pairs1[i].inten_sum;
	
	statistics.intensityInCharge2PairsIso2=0;
	statistics.numCharge2PairsIso2 = c2_pairs2.size();
	sort(c2_pairs2.begin(), c2_pairs2.end(),cmp_offset_pair_inten);
	maxIdx = (c2_pairs2.size()<4 ? c2_pairs2.size() : 4);
	for (size_t i=0; i<maxIdx; i++)
		statistics.intensityInCharge2PairsIso2+=c2_pairs2[i].inten_sum;
}







void PMCRankStats::clear()
{
	m_over_z=0;

	rank_score = NEG_INF;

	numFragmentPairs=0;
	numStrongFragmentPairs=0;
	numCharge2FragmentPairs=0;
	numStrongCharge2FragmentPairs=0;
	numH2oLossFragmentPairs=0;

	intensityInFragmentPairs=0;
	intensityInStrongPairs=0;
	intensityInCharge2Pairs=0;
	intensityInStrongCharge2Pairs=0;
	intensityInH2oLossPairs=0;
	intensityInCharge2H2oLossPairs=0;

	meanOffsetPairs=0;
	meanOffsetStrongPairs=0;
	meanOffsetCharge2Pairs=0;
	meanOffsetCharge2StrongPairs=0;
	meanOffsetH2oPairs=0;
	meanOffsetCharge2H2oPairs=0;

	ind_pairs_with_min_tol=false;			 
	ind_strong_pairs_with_min_tol=false;
	ind_c2_pairs_with_min_tol=false;
	ind_c2_strong_pairs_with_min_tol=false;

	offsetPairsOderedByIntensity.clear();
	strong_offsetPairsOderedByIntensity.clear();
	charge2offsetPairsOderedByIntensity.clear();


	numPairsIso0=0; intensityInPairsIso0=0;
	numPairsIso1=0; intensityInPairsIso1=0;
	numPairsIso2=0; intensityInPairsIso2=0;
}


void PMCSQS_Scorer::fillRankboostPmcSamples(const PeakList& pl,
								 int charge,
								 vector<MlSample>& samples) const
{
	const int numSamples = currentSpectrumPmcTables_[charge].size();
	const int idxsToSkip = int((1.0/bin_increment)+0.00001);
	vector<int> idxOffsets(4);
	
	idxOffsets[0]= -2*idxsToSkip;
	idxOffsets[1]= -idxsToSkip;
	idxOffsets[2]= idxsToSkip;
	idxOffsets[3]= 2*idxsToSkip;

	if (samples.size() != numSamples)
		samples.resize(numSamples);

	for (size_t i=0; i<numSamples; i++)
	{
		const PMCRankStats& stats = currentSpectrumPmcTables_[charge][i];
		MlSample& sam = samples[i];

		assert(currentSpectrumTotalIntensity_>=0.0);
		const float intensityNormalization = 1.0/(currentSpectrumTotalIntensity_+1.0);
		const mass_t mzOffset = (stats.m_over_z - pl.getHeader()->getMOverZ());

		size_t rankFeatureIdx=0;

		sam.clear();
		sam.addPair(rankFeatureIdx++,mzOffset);

		if (stats.numFragmentPairs<=2)
		{
			sam.addPair(rankFeatureIdx, mzOffset);
		}
		else if (stats.numFragmentPairs<4)
		{
			sam.addPair(rankFeatureIdx+1, mzOffset);
		}
		else
			sam.addPair(rankFeatureIdx+2, mzOffset);

		rankFeatureIdx+=3;

		if (stats.numStrongFragmentPairs<3)
		{
			sam.addPair(rankFeatureIdx, mzOffset);
		}
		else
			sam.addPair(rankFeatureIdx+1, mzOffset);

		rankFeatureIdx+=2;

		if (charge>=2)
		{
			if (stats.numCharge2FragmentPairs<=2)
			{
				sam.addPair(rankFeatureIdx, mzOffset);
			}
			else if (stats.numCharge2FragmentPairs<4)
			{
				sam.addPair(rankFeatureIdx+1, mzOffset);
			}
			else
				sam.addPair(rankFeatureIdx+2, mzOffset);

			rankFeatureIdx+=3;

			if (stats.numStrongCharge2FragmentPairs<3)
			{
				sam.addPair(rankFeatureIdx, mzOffset);
			}
			else
				sam.addPair(rankFeatureIdx+1, mzOffset);

			rankFeatureIdx+=2;
		}
		else
			rankFeatureIdx+=5;
			
	/*	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=2");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS >5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS <4");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS >4");

		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=2");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS >5");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS <4");
		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS >4");*/

		sam.addPair(rankFeatureIdx++,stats.numFragmentPairs);
		sam.addPair(rankFeatureIdx++,stats.numStrongFragmentPairs);

		if (charge>1)
		{
			sam.addPair(rankFeatureIdx++,stats.numCharge2FragmentPairs);
			sam.addPair(rankFeatureIdx++,stats.numStrongCharge2FragmentPairs);
		}
		else
			rankFeatureIdx+=2;

		sam.addPair(rankFeatureIdx++,stats.numH2oLossFragmentPairs);

		if (charge>1)
			sam.addPair(rankFeatureIdx,stats.numCharge2H2oLossFragmentPairs);
		rankFeatureIdx++;

		sam.addPair(rankFeatureIdx++,stats.intensityInFragmentPairs * intensityNormalization);
		sam.addPair(rankFeatureIdx++,stats.intensityInStrongPairs   * intensityNormalization);

		if (charge>1)
		{
			sam.addPair(rankFeatureIdx++,stats.intensityInCharge2Pairs  * intensityNormalization);
			sam.addPair(rankFeatureIdx++,stats.intensityInStrongCharge2Pairs  * intensityNormalization);
		}
		else
			rankFeatureIdx+=2;

		sam.addPair(rankFeatureIdx++,stats.intensityInH2oLossPairs        * intensityNormalization);

		if (charge>1)
			sam.addPair(rankFeatureIdx,stats.intensityInCharge2H2oLossPairs * intensityNormalization);
		rankFeatureIdx++;

		// averages of top k offsets
		float avg=0;
		for (size_t j =0; j<7 && j<stats.offsetPairsOderedByIntensity.size(); j++)
		{
			avg += fabs(stats.offsetPairsOderedByIntensity[j]);
			if (j>=2)
				sam.addPair(rankFeatureIdx+j-2,avg/static_cast<float>(j));
		}
		rankFeatureIdx+=5;

		if (charge>1)
		{
			avg=0;
			for (size_t j =0; j<7 && j<stats.charge2offsetPairsOderedByIntensity.size(); j++)
			{
				avg += fabs(stats.charge2offsetPairsOderedByIntensity[j]);
				if (j>=2)
					sam.addPair(rankFeatureIdx+j-2,avg/static_cast<float>(j));
			}
		}
		rankFeatureIdx+=5;


		// offset data
		if (stats.meanOffsetPairs<POS_INF)
		{
			sam.addPair(rankFeatureIdx++,stats.meanOffsetPairs);
			sam.addPair(rankFeatureIdx++,stats.meanOffsetPairs/(1.0+stats.numFragmentPairs));
		}
		else
			rankFeatureIdx+=2;

		if (stats.meanOffsetStrongPairs<POS_INF)
		{
			sam.addPair(rankFeatureIdx++,stats.meanOffsetStrongPairs);
			sam.addPair(rankFeatureIdx++,stats.meanOffsetStrongPairs/(1.0+stats.numStrongFragmentPairs));
		}
		else
			rankFeatureIdx+=2;

		if (charge>1 && stats.meanOffsetCharge2Pairs<POS_INF)
		{
			sam.addPair(rankFeatureIdx++,stats.meanOffsetCharge2Pairs);
			sam.addPair(rankFeatureIdx++,stats.meanOffsetCharge2Pairs/(1.0+stats.numCharge2FragmentPairs));
		}
		else
			rankFeatureIdx+=2;

		if (charge>1 && stats.meanOffsetCharge2StrongPairs<POS_INF)
		{
			sam.addPair(rankFeatureIdx++,stats.meanOffsetCharge2StrongPairs);
			sam.addPair(rankFeatureIdx++,stats.meanOffsetCharge2StrongPairs/(1.0+stats.numStrongCharge2FragmentPairs));
		}
		else
			rankFeatureIdx+=2;

		if (stats.meanOffsetH2oPairs<POS_INF)
		{
			sam.addPair(rankFeatureIdx++,stats.meanOffsetH2oPairs);
			sam.addPair(rankFeatureIdx++,stats.meanOffsetH2oPairs/(1.0+stats.numH2oLossFragmentPairs));
		}
		else
			rankFeatureIdx+=2;

		if (charge>1 && stats.meanOffsetCharge2H2oPairs<POS_INF)
		{
			sam.addPair(rankFeatureIdx++,stats.meanOffsetCharge2H2oPairs);
			sam.addPair(rankFeatureIdx++,stats.meanOffsetCharge2H2oPairs/(1.0+stats.numCharge2H2oLossFragmentPairs));
		}
		else
			rankFeatureIdx+=2;

		// individual offsets
		for (size_t j=0; j<5 && j<stats.offsetPairsOderedByIntensity.size(); j++)
			sam.addPair(rankFeatureIdx+j,stats.offsetPairsOderedByIntensity[j]);
		rankFeatureIdx+=5;

		if (charge>1)
			for (size_t j=0; j<5 && j<stats.charge2offsetPairsOderedByIntensity.size(); j++)
				sam.addPair(rankFeatureIdx+j,stats.charge2offsetPairsOderedByIntensity[j]);
		rankFeatureIdx+=5;

	
		// add the +0 +1 +2 strict counts
		sam.addPair(rankFeatureIdx++,stats.numPairsIso0);
		sam.addPair(rankFeatureIdx++,stats.intensityInPairsIso0 * intensityNormalization);
		sam.addPair(rankFeatureIdx++,stats.numPairsIso1);
		sam.addPair(rankFeatureIdx++,stats.intensityInPairsIso1 * intensityNormalization);
		sam.addPair(rankFeatureIdx++,stats.numPairsIso2);
		sam.addPair(rankFeatureIdx++,stats.intensityInPairsIso2 * intensityNormalization);
	
		if (charge>1)
		{
			sam.addPair(rankFeatureIdx++,stats.numCharge2PairsIso0);
			sam.addPair(rankFeatureIdx++,stats.intensityInCharge2PairsIso0 * intensityNormalization);
			sam.addPair(rankFeatureIdx++,stats.numCharge2PairsIso1);
			sam.addPair(rankFeatureIdx++,stats.intensityInCharge2PairsIso1 * intensityNormalization);
			sam.addPair(rankFeatureIdx++,stats.numCharge2PairsIso2);
			sam.addPair(rankFeatureIdx++,stats.intensityInCharge2PairsIso2 * intensityNormalization);
		}
		else
			rankFeatureIdx+=6;

		// add comparative features to -2 -1 +1 +2 Da away
		for (size_t j=0; j<idxOffsets.size(); j++)
		{
			const int other_idx = i + idxOffsets[j];
			if (other_idx<0 || other_idx>= samples.size())
			{
				rankFeatureIdx+=12;
				continue;
			}

			const PMCRankStats& other = currentSpectrumPmcTables_[charge][other_idx];

			sam.addPair(rankFeatureIdx++,stats.numFragmentPairs - other.numFragmentPairs);
			sam.addPair(rankFeatureIdx++,stats.numStrongFragmentPairs - other.numStrongFragmentPairs);

			if (charge>1)
			{
				sam.addPair(rankFeatureIdx++,stats.numCharge2FragmentPairs - other.numCharge2FragmentPairs);
				sam.addPair(rankFeatureIdx++,stats.numStrongCharge2FragmentPairs - other.numStrongCharge2FragmentPairs);
			}
			else
				rankFeatureIdx+=2;

			sam.addPair(rankFeatureIdx++,stats.numH2oLossFragmentPairs - other.numH2oLossFragmentPairs);

			if (charge>1)
				sam.addPair(rankFeatureIdx,stats.numCharge2H2oLossFragmentPairs - other.numCharge2H2oLossFragmentPairs);
			rankFeatureIdx++;

			sam.addPair(rankFeatureIdx++,(stats.intensityInFragmentPairs - other.intensityInFragmentPairs) * intensityNormalization);
			sam.addPair(rankFeatureIdx++,(stats.intensityInStrongPairs - other.intensityInStrongPairs) * intensityNormalization);
			sam.addPair(rankFeatureIdx++,(stats.intensityInCharge2Pairs - other.intensityInCharge2Pairs) * intensityNormalization);
			sam.addPair(rankFeatureIdx++,(stats.intensityInStrongCharge2Pairs - other.intensityInStrongCharge2Pairs) * intensityNormalization);
			sam.addPair(rankFeatureIdx++,(stats.intensityInH2oLossPairs - other.intensityInH2oLossPairs) * intensityNormalization);
			sam.addPair(rankFeatureIdx++,(stats.intensityInCharge2H2oLossPairs - other.intensityInCharge2H2oLossPairs) * intensityNormalization);
		}

		const int plus_idx = i + idxsToSkip;
		const int minus_idx = i-idxsToSkip;

		if (plus_idx<samples.size() && minus_idx>0)
		{
			const PMCRankStats& plus = currentSpectrumPmcTables_[charge][plus_idx];
			const PMCRankStats& minus = currentSpectrumPmcTables_[charge][minus_idx];

			sam.addPair(rankFeatureIdx++,plus.numFragmentPairs - minus.numFragmentPairs);
			sam.addPair(rankFeatureIdx++,plus.numStrongFragmentPairs - minus.numStrongFragmentPairs);
			sam.addPair(rankFeatureIdx++,plus.numCharge2FragmentPairs - minus.numCharge2FragmentPairs);
			sam.addPair(rankFeatureIdx++,plus.numStrongCharge2FragmentPairs - minus.numStrongCharge2FragmentPairs);
			sam.addPair(rankFeatureIdx++,plus.numH2oLossFragmentPairs - minus.numH2oLossFragmentPairs);
			sam.addPair(rankFeatureIdx++,plus.numCharge2H2oLossFragmentPairs - minus.numCharge2H2oLossFragmentPairs);

			sam.addPair(rankFeatureIdx++,(plus.intensityInFragmentPairs - minus.intensityInFragmentPairs) * intensityNormalization);
			sam.addPair(rankFeatureIdx++,(plus.intensityInStrongPairs - minus.intensityInStrongPairs) * intensityNormalization);
			sam.addPair(rankFeatureIdx++,(plus.intensityInCharge2Pairs - minus.intensityInCharge2Pairs) * intensityNormalization);
			sam.addPair(rankFeatureIdx++,(plus.intensityInStrongCharge2Pairs - minus.intensityInStrongCharge2Pairs) * intensityNormalization);
			sam.addPair(rankFeatureIdx++,(plus.intensityInH2oLossPairs - minus.intensityInH2oLossPairs) * intensityNormalization);
			sam.addPair(rankFeatureIdx++,(plus.intensityInCharge2H2oLossPairs - minus.intensityInCharge2H2oLossPairs) * intensityNormalization);
		}
	}
}




void init_PMC_feature_names(vector<string>& names)
{
	names.clear();
	int i;

	names.push_back("OFFSET FROM MEASURED M/Z");

	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=2");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=4");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS >4");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS <3");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS >=3");

	names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=2");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=4");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS >4");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS <3");
	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS >=3");
	
	names.push_back("# PAIRS");
	names.push_back("# STRONG PAIRS");
	names.push_back("# C2 PAIRS");
	names.push_back("# STRONG C2 PAIRS");
	names.push_back("# H2O PAIRS");
	names.push_back("# C2 H2O PAIRS");

	names.push_back("INTEN PAIRS");
	names.push_back("INTEN STRONG PAIRS");
	names.push_back("INTEN C2 PAIRS");
	names.push_back("INTEN STRONG C2 PAIRS");
	names.push_back("INTEN H2O PAIRS");
	names.push_back("INTEN C2 H2O PAIRS");

	for (i=2; i<7; i++)
	{
		char name[64];
		sprintf(name,"AVG OFFSET TOP (STRONG %d)",i);
		names.push_back(name);
	}

	for (i=2; i<7; i++)
	{
		char name[64];
		sprintf(name,"AVG OFFSET TOP C2 (STRONG %d)",i);
		names.push_back(name);
	}

	names.push_back("MEAN OFFSET PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET PAIRS");

	names.push_back("MEAN OFFSET STRONG PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET STRONG PAIRS");

	names.push_back("MEAN OFFSET C2 PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET C2 PAIRS");

	names.push_back("MEAN OFFSET STRONG C2 PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET STRONG C2 PAIRS");

	names.push_back("MEAN OFFSET H2O PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET H2O PAIRS");

	names.push_back("MEAN OFFSET C2 H2O PAIRS");
	names.push_back("WEIGHTED MEAN OFFSET C2 H2O PAIRS");

	for (i=0; i<5; i++)
	{
		char name[64];
		sprintf(name,"PAIR OFFSET (STRONG %d)",i+1);
		names.push_back(name);
	}

	for (i=0; i<5; i++)
	{
		char name[64];
		sprintf(name,"C2 PAIR OFFSET (STRONG %d)",i+1);
		names.push_back(name);
	}



	names.push_back("NUM STRICT 0");
	names.push_back("INTEN STRICT 0");

	names.push_back("NUM STRICT 1");
	names.push_back("INTEN STRICT 1");

	names.push_back("NUM STRICT 2");
	names.push_back("INTEN STRICT 2");

	names.push_back("NUM C2 STRICT 0");
	names.push_back("INTEN C2 STRICT 0");

	names.push_back("NUM C2 STRICT 1");
	names.push_back("INTEN C2 STRICT 1");

	names.push_back("NUM C2 STRICT 2");
	names.push_back("INTEN C2 STRICT 2");

	// diff features with -2 -1 +1 +2
	const string dis_labels[]={"-2","-1","+1","+2"};
	for (i=0; i<4; i++)
	{
		const string prefix = "DIFF WITH "+dis_labels[i]+" ";

		names.push_back(prefix+"# PAIRS");
		names.push_back(prefix+"# STRONG PAIRS");
		names.push_back(prefix+"# C2 PAIRS");
		names.push_back(prefix+"# STRONG C2 PAIRS");
		names.push_back(prefix+"# H2O PAIRS");
		names.push_back(prefix+"# C2 H2O PAIRS");

		names.push_back(prefix+"INTEN PAIRS");
		names.push_back(prefix+"INTEN STRONG PAIRS");
		names.push_back(prefix+"INTEN C2 PAIRS");
		names.push_back(prefix+"INTEN STRONG C2 PAIRS");
		names.push_back(prefix+"INTEN H2O PAIRS");
		names.push_back(prefix+"INTEN C2 H2O PAIRS");
	}

	names.push_back("DIFF +1/-1 # PAIRS");
	names.push_back("DIFF +1/-1 # STRONG PAIRS");
	names.push_back("DIFF +1/-1 # C2 PAIRS");
	names.push_back("DIFF +1/-1 # STRONG C2 PAIRS");
	names.push_back("DIFF +1/-1 # H2O PAIRS");
	names.push_back("DIFF +1/-1 # C2 H2O PAIRS");

	names.push_back("DIFF +1/-1 INTEN PAIRS");
	names.push_back("DIFF +1/-1 INTEN STRONG PAIRS");
	names.push_back("DIFF +1/-1 INTEN C2 PAIRS");
	names.push_back("DIFF +1/-1 INTEN STRONG C2 PAIRS");
	names.push_back("DIFF +1/-1 INTEN H2O PAIRS");
	names.push_back("DIFF +1/-1 INTEN C2 H2O PAIRS");
	cout << "Initialized: " << names.size() << " real feature names..." << endl;

	for (size_t i=0; i<names.size(); i++)
		cout << i << "\t" << names[i] << endl;
	exit(0);
}



/*
void PMCSQS_Scorer::output_pmc_rank_results(const FileManager& fm,
											int charge,
											const vector<SingleSpectrumFile *>& testHeaders) 
{
	BasicSpecReader bsr;
	static QCPeak peaks[5000];

	vector<int> org_offset_counts, new_offset_counts;
	org_offset_counts.resize(201,0);
	new_offset_counts.resize(201,0);

	vector<mass_t> org_offsets;
	vector<mass_t> corr_offsets;

	org_offsets.clear();
	corr_offsets.clear();

	int i;
	for (i=0; i<testHeaders.size(); i++)
	{
		SingleSpectrumFile* ssf = testHeaders[i];
		BasicSpectrum bs;
	
	//	bs.num_peaks = bsr.read_basic_spec(config_,fm,ssf,peaks);
		bs.peaks = peaks;
		bs.ssf = ssf;

	//	init_for_current_spec(config_,bs);
		calculate_curr_spec_pmc_values(bs, bin_increment);

		PmcSqsChargeRes res;
	//	find_best_mz_values_from_rank_model(bs, charge, config_->get_pm_tolerance(),res);

		ssf->peptide.calc_mass(config_);
		mass_t trueMz = (ssf->peptide.get_mass() + 18.01 + charge)/charge;

		org_offsets.push_back(trueMz - ssf->m_over_z);
		corr_offsets.push_back(trueMz - res.mz1);
	}

	mass_t m_org,sd_org,m_corr,sd_corr;
	calc_mean_sd(org_offsets,&m_org, &sd_org);
	calc_mean_sd(corr_offsets,&m_corr,&sd_corr);

	cout << "CHARGE: " << charge << endl;
	cout << "ORG:  mean " << m_org << " " << sd_org << endl;
	cout << "CORR: mean " << m_corr << " " << sd_corr << endl;

	for (i=0; i<org_offsets.size(); i++)
	{
		int org_idx = 100 + int(org_offsets[i] * 20);
		if (org_idx<0)
			org_idx = 0;
		if (org_idx>200)
			org_idx=200;
		org_offset_counts[org_idx]++;

		int new_idx = 100 + int(corr_offsets[i] * 20);
		if (new_idx<0)
			new_idx = 0;
		if (new_idx>200)
			new_idx=200;
		new_offset_counts[new_idx]++;
	}

	int cum_org=0;
	int cum_new=0;
	for (i=0; i<=200; i++)
	{

		if (org_offset_counts[i]==0 && new_offset_counts[i]==0)
			continue;
		
		cum_org+=org_offset_counts[i];
		cum_new+=new_offset_counts[i];
		cout << fixed << setprecision(3) << i*0.05 - 5.0 << "\t" <<
			org_offset_counts[i]/(float)org_offsets.size() << "\t" <<
			new_offset_counts[i]/(float)corr_offsets.size() << "\t" <<
			cum_org/(float)org_offsets.size() << "\t"<<
			cum_new/(float)corr_offsets.size() << endl;

	}


}
*/

/**********************************************************************************
Finds the best m/z values for a given charge
If the pm_tolerance is low (less than 0.5) then the m/z vlaues need to reflect +-X Da
from the recorded m/z
***********************************************************************************/
void PMCSQS_Scorer::computeBestMzValuesForCharge(const PeakList& pl, 
	     						  int charge,
								  mass_t precursorMassTolerance,
								  PmcSqsChargeRes& result)
{
		if (charge<1)
	{
		cout << "Error: trying to find m/z values which carge < 1 !" << endl;
		exit(1);
	}

	static vector<RankBoostSample> spectrumRankSamples;
	static vector<float> rankingScores;

	const mass_t spectrumMOverZ  = pl.getHeader()->getMOverZ();
	const mass_t maxMzDiffAllowed = precursorMassTolerance / static_cast<mass_t>(charge);

	fillRankboostPmcSamplesOld(pl, charge, spectrumRankSamples);
//	fill_RankBoost_smaples_with_PMC(bs, charge, spectrumRankSamples);

	if (rankingScores.size() != spectrumRankSamples.size())
		rankingScores.resize(spectrumRankSamples.size(),NEG_INF);

	const mass_t pm_with_19 = spectrumMOverZ * charge - (charge + 1);
	const size_t sizeIndex = calcRankModelSizeIndex(charge, pm_with_19);
	
	int best_idx=-1;
	float best_score=NEG_INF;

	if (charge>= pmc_rank_models.size() ||
		sizeIndex>= pmc_rank_models[charge].size() ||
		! pmc_rank_models[charge][sizeIndex])
	{
		//
	}
	else
	{
	//	cout << "spec: " << spectrumMOverZ << "  charge: " << charge << endl;
		int i;
		for (i=0; i<spectrumRankSamples.size(); i++)
		{
			// make sure the m/z is in an allowed range (i.e., it has a mass that can possibly be a +-X Da
			// shift from the true pm). Assume shift can be at most
			if (precursorMassTolerance<0.5)
			{
				const mass_t table_mz = currentSpectrumPmcTables_[charge][i].m_over_z;
				const mass_t one_over_charge = 1.0033 / (mass_t)charge; // ~mass of isotopic peak difference
				int d;
				for (d=-4; d<=4; d++)
				{
					const mass_t mz_diff = fabs(table_mz + d*one_over_charge - spectrumMOverZ);
					if (mz_diff<maxMzDiffAllowed)
						break;
				}

				if (d>4)
				{
					rankingScores[i] = NEG_INF;
					continue;
				}
			//	cout << "ok : " << table_mz << endl;
			}

			rankingScores[i]=pmc_rank_models[charge][sizeIndex]->calc_rank_score(spectrumRankSamples[i]);
			currentSpectrumPmcTables_[charge][i].rank_score = rankingScores[i];
			if (rankingScores[i]>best_score)
			{
				best_score=rankingScores[i];
				best_idx = i;
			}
		}
	}

	// no suitable models were found for this spectrum
	if (best_idx<0)
	{
		result.mz1 = pl.getHeader()->getMOverZ();
		result.score1 = 10.0;
		result.mz2 = NEG_INF;
		result.score2 = NEG_INF;
		return;
	}

	
	result.mz1 = currentSpectrumPmcTables_[charge][best_idx].m_over_z;
	result.score1 = best_score;

	// look for additional m/z
	int second_best_idx=-1;
	float second_best_score=NEG_INF;

	const mass_t mz_diff = currentSpectrumPmcTables_[charge][1].m_over_z - 
						   currentSpectrumPmcTables_[charge][0].m_over_z;

	const int idx_diff = (int)(0.45/(charge * mz_diff));

	int i;
	for (i=0; i<spectrumRankSamples.size(); i++)
	{
		if (rankingScores[i]>NEG_INF && abs(i-best_idx)<idx_diff)
			continue;

		if (rankingScores[i]>second_best_score)
		{
			second_best_score=rankingScores[i];
			second_best_idx = i;
		}

	}
 
	if (second_best_idx>=0)
	{
		result.mz2 = currentSpectrumPmcTables_[charge][second_best_idx].m_over_z;
		result.score2 = second_best_score;
	} 
	else
	{
		result.mz2 = NEG_INF;
		result.score2 = NEG_INF;
	}
}

/**********************************************************************************
Finds the best m/z values for a given charge
If the pm_tolerance is low (less than 0.5) then the m/z vlaues need to reflect +-X Da
from the recorded m/z
***********************************************************************************
void PMCSQS_Scorer::find_best_mz_values_from_rank_model(
										const BasicSpectrum& bs, 
										int charge,
										mass_t pm_tolerance,
										PmcSqsChargeRes& res)
{
	static vector<RankBoostSample> spectrumRankSamples;
	static vector<float> rankingScores;

	const mass_t spectrumMOverZ = bs.ssf->m_over_z;
	const mass_t maxMzDiffAllowed = pm_tolerance / charge;

	fill_RankBoost_smaples_with_PMC(bs, charge, spectrumRankSamples);

	if (rankingScores.size() != spectrumRankSamples.size())
		rankingScores.resize(spectrumRankSamples.size(),NEG_INF);

	const mass_t pm_with_19 = bs.ssf->m_over_z * charge - (charge + 1);
	const size_t sizeIndex = calcRankModelSizeIndex(charge, pm_with_19);
	
	int best_idx=-1;
	float best_score=NEG_INF;

	if (charge>= pmc_rank_models.size() ||
		sizeIndex>= pmc_rank_models[charge].size() ||
		! pmc_rank_models[charge][sizeIndex])
	{
		//
	}
	else
	{
	//	cout << "spec: " << spectrumMOverZ << "  charge: " << charge << endl;
		int i;
		for (i=0; i<spectrumRankSamples.size(); i++)
		{
			// make sure the m/z is in an allowed range (i.e., it has a mass that can possibly be a +-X Da
			// shift from the true pm). Assume shift can be at most
			if (pm_tolerance<0.5)
			{
				const mass_t table_mz = currentSpectrumPmcTables_[charge][i].m_over_z;
				const mass_t one_over_charge = 1.0033 / (mass_t)charge; // ~mass of isotopic peak difference
				int d;
				for (d=-4; d<=4; d++)
				{
					const mass_t mz_diff = fabs(table_mz + d*one_over_charge - spectrumMOverZ);
					if (mz_diff<maxMzDiffAllowed)
						break;
				}

				if (d>4)
				{
					rankingScores[i] = NEG_INF;
					continue;
				}
			//	cout << "ok : " << table_mz << endl;
			}

			rankingScores[i]=pmc_rank_models[charge][sizeIndex]->calc_rank_score(spectrumRankSamples[i]);
			currentSpectrumPmcTables_[charge][i].rank_score = rankingScores[i];
			if (rankingScores[i]>best_score)
			{
				best_score=rankingScores[i];
				best_idx = i;
			}
		}
	}

	// no suitable models were found for this spectrum
	if (best_idx<0)
	{
		res.mz1 = bs.ssf->m_over_z;
		res.score1 = 10.0;
		res.mz2 = NEG_INF;
		res.score2 = NEG_INF;
		return;
	}

	
	res.mz1 = currentSpectrumPmcTables_[charge][best_idx].m_over_z;
	res.score1 = best_score;

	// look for additional m/z
	int second_best_idx=-1;
	float second_best_score=NEG_INF;

	const mass_t mz_diff = currentSpectrumPmcTables_[charge][1].m_over_z - 
						   currentSpectrumPmcTables_[charge][0].m_over_z;

	const int idx_diff = (int)(0.45/(charge * mz_diff));

	int i;
	for (i=0; i<spectrumRankSamples.size(); i++)
	{
		if (rankingScores[i]>NEG_INF && abs(i-best_idx)<idx_diff)
			continue;

		if (rankingScores[i]>second_best_score)
		{
			second_best_score=rankingScores[i];
			second_best_idx = i;
		}

	}
 
	if (second_best_idx>=0)
	{
		res.mz2 = currentSpectrumPmcTables_[charge][second_best_idx].m_over_z;
		res.score2 = second_best_score;
	} 
	else
	{
		res.mz2 = NEG_INF;
		res.score2 = NEG_INF;
	}

//	const int sizeIndex = this->get_rank_model_sizeIndex(charge,res.mz1*charge-charge+1);
//	res.mz1 -= pmcMzBiases_[charge][sizeIndex];
//	res.mz2 -= pmcMzBiases_[charge][sizeIndex];

//	cout << charge << " ]\t" << res.mz1 << "\t" << res.prob1 << "\t" << res.mz2 << "\t" << res.prob2 << endl;


}*/


void PMCSQS_Scorer::benchmarkPm(const Config* config, const char* list)
{
	SpectraAggregator sa;
	sa.initializeFromTextFile(list, config);

	SpectraList sl(sa);
	sl.selectAllAggregatorHeaders();
	
	cout << "Found " << sl.getNumHeaders() << " spectra...";

	if (sl.getNumHeaders() == 0)
		return;

	int numFilesWithoutSpectra = 0;
	size_t peakBufferSize = 10000;
	Peak*  peakBuffer = new Peak[peakBufferSize];

	const vector< vector<mass_t> >& sizes = this->pmcMassThresholds_;
	
	vector< vector<int> > numTested(sizes.size()); // charge / size
	vector< vector<int> > numWrongCharge(sizes.size());
	vector< vector< vector<mass_t> > > deltas(sizes.size());
	int totalTested=0;

	for (size_t i=0; i<sizes.size(); i++)
	{
		const size_t numSizes = sizes[i].size()+1;
		numTested[i].resize(numSizes,0);
		numWrongCharge[i].resize(numSizes,0);
		deltas[i].resize(numSizes);
	}

	for (size_t j=0; j<sl.getNumHeaders(); j++)
	{
		const SingleSpectrumHeader* header = sl.getSpectrumHeader(j);
		if (header->getOriginalNumPeaks()>1e6)
			continue;

		if (header->getPeptideStr().length()<3)
			continue;

		if (header->getOriginalNumPeaks()>= peakBufferSize)
		{
			delete [] peakBuffer;
			peakBufferSize = header->getOriginalNumPeaks()*2;
			peakBuffer = new Peak[peakBufferSize];
		}

		PeakList pl;
		pl.setPeaksPtr( peakBuffer );

        // NP3 GOT change few peaks from 7 to 1 @@@
		if (pl.readPeaksToBuffer(sa, header, peakBuffer) < 1) // if not enough peaks read, skip this spectrum
			continue;

		pl.initializePeakList(config_, true);

		// NP3 GOT change few peaks from 7 to 1 @@@
		if (pl.getNumPeaks()<1) // don't bother with spectra with too few peaks
			continue;

		size_t maxCharge=0;
		const float sqs = calculateSqsScore(config_, pl, &maxCharge);

		PmcSqsChargeRes res;
		computeBestMzValuesForCharge(pl, maxCharge, config_->get_pm_tolerance(), res);

		Peptide pep;
		pep.parseFromString(config, header->getPeptideStr());
		const int trueCharge = static_cast<int>(header->getCharge());
		const size_t sizeIdx = calcRankModelSizeIndex(trueCharge, pep.get_mass_with_19());

		numTested[trueCharge][sizeIdx]++;
		totalTested++;
		if (trueCharge != maxCharge)
		{
			numWrongCharge[trueCharge][sizeIdx]++;
			continue;
		}

		const mass_t trueMz =  (pep.get_mass_with_19() + (trueCharge-1.0)*MASS_PROTON) / trueCharge;
		mass_t delta = trueMz - res.mz1;
		deltas[trueCharge][sizeIdx].push_back(delta);
	//	cout << header->getMOverZ() << " : " << maxCharge << "\t" << res.mz1 << "\t" << delta << endl;
	}

	cout << endl << "Summary:" << endl;
	cout << "--------" << endl;
	cout << fixed << setprecision(3);
	cout << "Benchmarked against " << totalTested << " spectra." << endl;
	for (size_t i=0; i<sizes.size(); i++)
	{
		for (size_t j=0; j<=sizes[i].size(); j++)
		{
			if (numTested[i][j]<1)
				continue;
			cout << "charge " << i << "\tsize " << j << "\t#tested\t" << numTested[i][j] << "\t";
			cout << " wrong " << numWrongCharge[i][j] << " (" 
				 << static_cast<float>(numWrongCharge[i][j])/static_cast<float>(numTested[i][j]) <<")\t";

			if (deltas[i][j].size()<3)
				continue;
		
			MeanSdStats msd;
			int numWithin02=0;
			int numWithin03=0;
			int numWithin04=0;
			int numWithin05=0;

			for (size_t k=0; k<deltas[i][j].size(); k++)
			{
				mass_t d=deltas[i][j][k];
				if (d<0.2)
					numWithin02++;
				if (d<0.3)
					numWithin03++;
				if (d<0.4)
					numWithin04++;
				if (d<0.5)
					numWithin05++;
				msd.addWX(1.0,d);
			}

			float n = static_cast<float>(deltas[i][j].size());
			double m,sd;
			msd.calcMeanAndSd(m,sd);
			cout << "M: " << m << "\tSD: " << sd << "\t" << "0.2: " << numWithin02/n << " 0.3: "
				<< numWithin03/n << " 0.4:" << numWithin04/n << " 0.5:" << numWithin05/n << endl;
		}
	}
}






