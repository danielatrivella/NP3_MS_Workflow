#include "PMCSQS.h"
#include "SpectraList.h"
#include "PepNovo_auxfun.h"
#include <assert.h>


void PMCSQS_Scorer::setClassWeightsAccordingToData(const SpectraAggregator& positiveSpectra,
												   const SpectraAggregator& negativeSpectra,
												   double	weightNegativeSpectra,
												   vector< vector<double> >& classWeights,
												   bool verbose) const
{
	const size_t numSizes = static_cast<size_t>(getNumSizes());
	vector< vector<int> > counters(numSizes, vector<int>(maximalChargeWithModels_ + 1,0));// size / charge

	classWeights.clear();
	classWeights.resize(numSizes, vector<double>(maximalChargeWithModels_ + 1,0.0));

	assert(weightNegativeSpectra>0.01 && weightNegativeSpectra<0.99);

	// count weights of negative spectra
	SpectraList slNegative(negativeSpectra);
	slNegative.selectAllAggregatorHeaders();
	for (size_t i=0; i<slNegative.getNumHeaders(); i++)
	{
		const int sizeIndex = getSqsSizeIndex(slNegative.getSpectrumHeader(i)->getMOverZ());
		counters[sizeIndex][0]++;
	}

	// count weights of positive spectra
	SpectraList slPositive(positiveSpectra);
	slPositive.selectAllAggregatorHeaders();
	for (size_t i=0; i<slPositive.getNumHeaders(); i++)
	{
		const SingleSpectrumHeader* header = slPositive.getSpectrumHeader(i);
		const int sizeIndex = getSqsSizeIndex(header->getMOverZ());
		const int charge	= header->getCharge();
		if (charge<=0 || charge>=counters[sizeIndex].size())
		{
			cout << "Warning: positive training spectrum has charge " << charge << endl;
			continue;
		}
		counters[sizeIndex][charge]++;
	}


	if (verbose)
	{
		cout << "Found following spectra for training SQS:" << endl;
		for (int i=0; i<counters.size(); i++)
		{
			cout << "size " << i;
			for (int j=0; j<counters[i].size(); j++)
				cout << "\t" << counters[i][j];
			cout << endl;
		}
		cout << endl;
	}

	// compute size distribution for negative spectra
	int totalNegativeCounts=0;
	for (size_t sizeIndex=0; sizeIndex<counters.size(); sizeIndex++)
		totalNegativeCounts+=counters[sizeIndex][0];

	for (size_t sizeIndex=0; sizeIndex<counters.size(); sizeIndex++)
		classWeights[sizeIndex][0] = weightNegativeSpectra * 
		(static_cast<double>(counters[sizeIndex][0])/static_cast<double>(totalNegativeCounts));

	int totalPositiveSpectra=0;
	for (size_t sizeIndex=0; sizeIndex<counters.size(); sizeIndex++)
		for (size_t charge=1; charge<counters[sizeIndex].size(); charge++)
			totalPositiveSpectra+= counters[sizeIndex][charge];

	for (size_t sizeIndex=0; sizeIndex<counters.size(); sizeIndex++)
		for (size_t charge=1; charge<counters[sizeIndex].size(); charge++)
			classWeights[sizeIndex][charge] = (1.0 - weightNegativeSpectra) *
				static_cast<double>(counters[sizeIndex][charge])/totalPositiveSpectra;

	if (verbose)
	{
		cout << "Assigned weights for SQS models as follows:" << endl;
		cout << "\t		Charges (0=noise)" << endl;
		cout << "Size";
		for (size_t c=0; c<classWeights[0].size(); c++)
			cout << "\t" << c;
		cout << endl;
		for (size_t sizeIndex=0; sizeIndex<classWeights.size(); sizeIndex++)
		{
			cout << sizeIndex << setprecision(4) << fixed;
			for (size_t c=0; c<classWeights[sizeIndex].size(); c++)
				cout << "\t" << classWeights[sizeIndex][c];
			cout << endl;
		}
	}
}





void PMCSQS_Scorer::trainRegressionSqsModels(const Config* config, 
											 const SpectraAggregator& positiveSpectra,
											 const char* pathNegativeSpectraList,
											 double weightNegativeSpectra)
{

	// initialize models and data structs

	const size_t maxHeadersPerCharge = 50000;
	maximalChargeWithModels_ = positiveSpectra.getMaxSepctrumCharge();

	setFragmentPairSumOffset(MASS_PROTON); // b+y - PM+19
	setBinMassIncrement(0.1);

	set_sqs_mass_thresholds();
	if (pmcMassThresholds_.size() == 0)
		pmcMassThresholds_=config->get_size_thresholds();

	const size_t numSizes = sqsMassThresholds_.size()+1;
	cout << "number of sizes for SQS models " << numSizes << endl;

	regressionSqsModels_.resize(numSizes,0);
	regressionChargeModels_.resize(numSizes);
	for (size_t i=0; i<numSizes; i++)
	{
		regressionChargeModels_[i].resize(maximalChargeWithModels_+1);
		for (size_t j=0; j<=maximalChargeWithModels_; j++)
			regressionChargeModels_[i][j].resize(maximalChargeWithModels_+1,0);
	}
	sqsRegressionNormalizingConstants_.resize(numSizes,0.0);

	SpectraAggregator negativeSpectra;
	negativeSpectra.initializeFromTextFile(pathNegativeSpectraList, config);

	vector< vector<double> > classWeights;
	setClassWeightsAccordingToData(positiveSpectra, negativeSpectra, weightNegativeSpectra, classWeights);
	string inputName = config->getModelName() + "_SQS";
	
	// create models according to size

	for (size_t sizeIndex=0; sizeIndex<numSizes; sizeIndex++)
	{	
		cout << endl << "Collecting training data for size " << sizeIndex << endl;

		const mass_t minMass = (sizeIndex == 0 ? 0 : sqsMassThresholds_[sizeIndex-1]);
		const mass_t maxMass = (sizeIndex == numSizes-1 ? MAX_FLOAT : sqsMassThresholds_[sizeIndex]);

		SpectraList negativeSl(negativeSpectra);
		SpectraList positiveSl(positiveSpectra);
		negativeSl.selectHeaders(minMass, maxMass);

		cout << "Found " << negativeSl.getNumHeaders() << " crap spectra in range " << minMass << " to " << maxMass << endl;
		
		for (size_t charge=0; charge<=maximalChargeWithModels_; charge++)
		{
			char trainingFilePath[256];
			sprintf(trainingFilePath,"%s/%s_SQS_%d_%d_dat.txt",config->getTrainingDir().c_str(),
															   config->getModelName().c_str(),
															   sizeIndex, charge);
			// put the file creation in a separate scope
			ifstream trainingStream(trainingFilePath);

			if (! trainingStream.good())
			{
				trainingStream.close();
				MlDataSet mld;
				mld.reserve(maxHeadersPerCharge);

				const SpectraAggregator& sa = (charge == 0 ? negativeSpectra : positiveSpectra);
				SpectraList sl(sa);

				if (charge>0)
				{
					sl.selectHeaders(minMass, maxMass, charge, charge);
				}
				else
					sl.selectHeaders(minMass, maxMass);

				if (sl.getNumHeaders()>maxHeadersPerCharge)
				{
					cout << "Reducing to " << maxHeadersPerCharge << endl;
					sl.randomlyReduceListToSize(maxHeadersPerCharge);
				}

				cout << "	charge " << charge << ", " << sl.getNumHeaders() << " spectra ";
					
				size_t tenth = 1 + (sl.getNumHeaders() / 10);
				for (size_t i=0; i<sl.getNumHeaders(); i++)
				{
					const SingleSpectrumHeader* header = sl.getSpectrumHeader(i);
					PeakList pl;
					MlSample sample;

					pl.readPeaksToLocalAllocation(sa,header);
					pl.initializePeakList(config, true);
					initializeForCurrentSpectrum(config, pl);
					calculatePmcValuesForSqs(pl, bin_increment);
					sample.label = charge;
					fillSqsSample(pl, sample);
					mld.addSample(sample);
					if (i % tenth == 0)
						cout << ".";
				}
				cout << endl;
				mld.printDatasetStatistics();	
				mld.writeWholeDataFile(trainingFilePath);
				cout << "Wrote " << trainingFilePath << endl;
			}
			else
			{
				cout << "Found " << trainingFilePath << endl;
				trainingStream.close();
			}
		}


		// model for all good spectra against bad
		if (1)
		{
			cout << endl;
			cout << "-----------------------------------------------------------------" << endl;
			cout << "Size " << sizeIndex << ": Training model for good spectra (all charges) vs. bad:" << endl;
			cout << "-----------------------------------------------------------------" << endl << endl;

			vector<double> weights(2,0.0);
			weights[0]=classWeights[sizeIndex][0];
			for (size_t charge=1; charge<=maximalChargeWithModels_; charge++)
				weights[1]+=classWeights[sizeIndex][charge];

			char classZeroFilePath[128];
			ostringstream oss;
			oss << "_" << sizeIndex << "_0";
			string outputName = inputName + oss.str();
			sprintf(classZeroFilePath,"%s/%s_SQS_%d_0_dat.txt",config->getTrainingDir().c_str(),
														   config->getModelName().c_str(),sizeIndex);
			MlTrainingContainer params;
			params.setTrainingAlgorithm(SMT_LOGISTIC_CG);
			params.addTrainingSetPath(classZeroFilePath);

			// add good spectra from all charges
			for (size_t charge=1; charge<=maximalChargeWithModels_; charge++)
			{
				char otherFilePath[128];
				sprintf(otherFilePath,"%s/%s_SQS_%d_%d_dat.txt",config->getTrainingDir().c_str(),
															config->getModelName().c_str(),sizeIndex, charge);
				params.addTrainingSetPath(otherFilePath);
			}
				
			params.setTestRatio(0.333);
			params.setInputDir(config->getTrainingDir().c_str());
			params.setInputName(inputName.c_str());
			params.setOutputName(outputName.c_str());
			params.setMaxNumIterations(200);
			params.setLambda(0.01);
			params.setVerboseLevel(1);
			params.setIndReadPreviousModel(true);
			params.setDesiredClassWeights(weights);

			char logFilePath[128];
			sprintf(logFilePath,"%s/%s_SQS_%d_0_log.txt",config->getTrainingDir().c_str(),
														 config->getModelName().c_str(),sizeIndex);

			regressionSqsModels_[sizeIndex] = new MlModel;
			MlModel *mlModel = regressionSqsModels_[sizeIndex];
			mlModel->createModel(params);
			mlModel->outputScoreModelAnalysisToLog(&params, logFilePath);

			// filter dataset to remove samples of class 0 that have a high probability of not being noise
			// (anything from class 0 that has less than 50% is likely to be an unidentified peptide that
			// got thrown in)

			vector<IdxVal> pairs;
			const vector<MlSample>& samples = params.getTrainingDataSet()->getSamples();
			const MlScoreModel* scoreModel = mlModel->getScoreModel();
			for (size_t i=0; i<samples.size(); i++)
			{
				if (samples[i].label>0)
						continue;
				double prob=scoreModel->calcScoreForClass(&samples[i],0);
				pairs.push_back(IdxVal(i,static_cast<value_t>(prob)));
			}
			sort(pairs.begin(),pairs.end());
			vector<size_t> idxsToRemove;
			for (size_t i=pairs.size()-1; i>0; i--)
				if (pairs[i].value<0.5)
					idxsToRemove.push_back(pairs[i].index);

			if (static_cast<double>(idxsToRemove.size())/static_cast<double>(pairs.size()) > 0.02)
			{
				cout << endl << "Removing " << idxsToRemove.size() << "/" << pairs.size() << " from class 0." << endl;

				params.removeTrainingAndTestIdxs(idxsToRemove);
				params.setIndReadPreviousModel(false);
				params.setIndClearWeights(false);
				params.setMaxNumIterations(300);

				cout << endl << "Retraining model... (" << params.getTrainingIdxs().size() 
					<< " , " << params.getTestIdxs().size() << ")" << endl << endl;

				mlModel->createModel(params);
				mlModel->outputScoreModelAnalysisToLog(&params, logFilePath);
			}

			// set normalization constant
			sqsRegressionNormalizingConstants_[sizeIndex] = 
				findSqsNormalizingConstant(sizeIndex, &params);

			cout << endl << "Performance using SQS model all good spectra against bad: " << endl;
			cout << endl << "Charge\tWeight\t%Bad \t%Good" << endl;
			cout		 << "------\t------\t-----\t-----" << endl;
			for (size_t charge=0; charge<=maximalChargeWithModels_; charge++)
			{
				char filePath[256];
				sprintf(filePath,"%s/%s_SQS_%d_%d_dat.txt",config->getTrainingDir().c_str(),
														  config->getModelName().c_str(),sizeIndex, charge);
				MlDataSet mld;
				mld.readDataFile(filePath);
				int numAssignedToBad=0;
				for (size_t i=0; i<mld.getNumSamples(); i++)
				{
					double p=calculateSqsScoreFromSample(mld.getSample(i), sizeIndex);
					if (p<0.1)
						numAssignedToBad++;
				}
				double badRatio = static_cast<double>(numAssignedToBad)/mld.getNumSamples();
				cout << charge << "\t" << classWeights[sizeIndex][charge] 
					 << "\t" << badRatio << "\t" << 1.0 - badRatio << endl;		
			}
			cout << endl;
		}

		// train models for charge1 vs. charge 2
		for (size_t charge1=1; charge1<maximalChargeWithModels_; charge1++)
		{
			for (size_t charge2=charge1+1; charge2<=maximalChargeWithModels_; charge2++)
			{
				char pathDat1[256], pathDat2[256], logFilePath[256];
				sprintf(pathDat1,"%s/%s_SQS_%d_%d_dat.txt",config->getTrainingDir().c_str(),
															config->getModelName().c_str(),sizeIndex, charge1);
				sprintf(pathDat2,"%s/%s_SQS_%d_%d_dat.txt",config->getTrainingDir().c_str(),
															config->getModelName().c_str(),sizeIndex, charge2);
				sprintf(logFilePath,"%s/%s_SQS_%d_%d_v_%d_log.txt",config->getTrainingDir().c_str(),
														  config->getModelName().c_str(),sizeIndex, charge1, charge2);
				vector<weight_t> weights(2);
				weights[0]=classWeights[sizeIndex][charge1];
				weights[1]=classWeights[sizeIndex][charge2];

				if (weights[0]<0.0005 || weights[1]< 0.0005)
					continue;
			
				ostringstream oss;
				oss << "_" << sizeIndex << "_" << charge1 << "_v_" << charge2;
				
				string outputName = inputName + oss.str();

				MlTrainingContainer params;
				params.setTrainingAlgorithm(SMT_LOGISTIC_CG);
				params.addTrainingSetPath(pathDat1);
				params.addTrainingSetPath(pathDat2);
				params.setTestRatio(0.333);
				params.setInputDir(config->getTrainingDir().c_str());
				params.setInputName(inputName.c_str());
				params.setOutputName(outputName.c_str());
				params.setMaxNumIterations(200);
				params.setLambda(0.01);
				params.setVerboseLevel(1);
				params.setIndReadPreviousModel(true);
				params.setDesiredClassWeights(weights);
				params.setLabelToConvertTo0(charge1);
			
				cout << endl << "*************************************************************" << endl;
				cout << endl << "Training model for size " << sizeIndex << " charge " << charge1 << " vs. charge " << charge2 << endl;
				cout << endl << "*************************************************************" << endl << endl;

				regressionChargeModels_[sizeIndex][charge1][charge2] = new MlModel;
				MlModel *mlModel = regressionChargeModels_[sizeIndex][charge1][charge2];
				mlModel->createModel(params);
				mlModel->outputScoreModelAnalysisToLog(&params, logFilePath);

				// test model
				const MlDataSet* testSet = params.getTestDataSet();
				const vector<weight_t>& testWeights = testSet->getClassWeights();
				vector < vector<weight_t> > assignments(2, vector<weight_t>(2,0.0));
				for (size_t i=0; i<testSet->getNumSamples(); i++)
				{
					const MlSample& sample = testSet->getSample(i);
					double p=mlModel->getScoreModel()->calcScoreForClass(&sample,0);
					if (p>0.5)
					{
						if (sample.label == 0) // label 0 == charge1
						{
							assignments[0][0]+=sample.weight;
						}
						else
							assignments[0][1]+=sample.weight;
					}
					else
					{
						if (sample.label == 0) // label 0 == charge1
						{
							assignments[1][0]+=sample.weight;
						}
						else
							assignments[1][1]+=sample.weight;
					}
				}
				const weight_t weightsClass0 = assignments[0][0]+assignments[0][1];
				const weight_t weightsClass1 = assignments[1][0]+assignments[1][1];
				const weight_t totalWeight = weightsClass0 + weightsClass1;

				cout << endl << "Charge assignment perofrmance:" << endl;
				cout <<         "------------------------------" << endl;
				cout << fixed << setprecision(3) << "\tWeight\tCh. " << charge1 << "\tCh. " << charge2 << endl;
				cout << "Ch. " << charge1 << "\t" << weightsClass0/totalWeight << "\t" 
					<< assignments[0][0]/weightsClass0 << "\t" << assignments[0][1]/weightsClass0 << endl;
				cout << "Ch. " << charge2 << "\t" << weightsClass1/totalWeight << "\t" 
					<< assignments[1][0]/weightsClass1 << "\t" << assignments[1][1]/weightsClass1 << endl;
				cout << endl << "Proportion correctly assigned: " 
					<< (assignments[0][0] + assignments[1][1])/totalWeight << endl;
				cout << endl;	
			}
		}
	}

	indInitializedSqs_ = true;

	string path;
}


struct SqsTuple {
	bool operator< (const SqsTuple& rhs) const
	{
		return (prob<rhs.prob);
	}

	int label;
	weight_t weight;
	double prob;
};

struct CumSqsTuple {
	float sqs;
	int numGoodBelow;
	int numBadBelow;
};



double PMCSQS_Scorer::findSqsNormalizingConstant(size_t sizeIndex,
												 MlTrainingContainer* params, 
												 double ratioPositiveMissed,
												 bool verbose)
{
	double threshold = 0.0;
	if (sizeIndex>=regressionSqsModels_.size() ||
		! regressionSqsModels_[sizeIndex]->getScoreModel() ||
		! regressionSqsModels_[sizeIndex]->getScoreModel()->getIndWasInitialized() )
		error("findSqsNormalizingConstant run on uninitialized model!");

	const MlDataSet* ds = params->getTestDataSet();
	const vector<size_t>& idxs = params->getTestIdxs();
	const MlScoreModel* scoreModel = regressionSqsModels_[sizeIndex]->getScoreModel();

	vector<SqsTuple> tuples(idxs.size());
	weight_t totalZero=0.0, totalRest = 0.0;
	for (size_t i=0; i<idxs.size(); i++)
	{
		const MlSample& sample = ds->getSample(idxs[i]);
		double p=scoreModel->calcScoreForClass(&sample,0);
		tuples[i].label = sample.label;
		tuples[i].weight= sample.weight;
		tuples[i].prob  = 1.0 - p;
		if (sample.label == 0)
		{
			totalZero += sample.weight;
		}
		else
			totalRest += sample.weight;

	//	cout << tuples[i].label << "\t" << tuples[i].prob << "\t" << tuples[i].weight << endl;
	}
	sort(tuples.begin(),tuples.end());

	const weight_t thresh = ratioPositiveMissed * totalRest;
	weight_t wZero=0.0, wRest=0.0;
	size_t i;
	for (i=0; i<tuples.size() && wRest<thresh; i++)
	{
		if (tuples[i].label == 0)
		{
			wZero += tuples[i].weight;
		}
		else
			wRest += tuples[i].weight;
	}
	threshold=tuples[i].prob;

	cout << "Thresh is: " << threshold << endl;
	cout << "Ratio positive above : " << 1.0 - wRest/totalRest << endl;
	cout << "Ratio negative above : " << 1.0 - wZero/totalZero << endl;

	// make histogram
	if (verbose)
	{
		vector<double> ratioZero(101,0.0);
		vector<double> ratioRest(101,0.0);
		weight_t wZero=0.0, wRest=0.0;
		for (size_t i=0; i<tuples.size(); i++)
		{
			if (tuples[i].label == 0)
			{
				wZero += tuples[i].weight;
			}
			else
				wRest += tuples[i].weight;
		
			const size_t idx = static_cast<size_t>(tuples[i].prob * 100.0);
			if (ratioZero[idx] == 0.0 || ratioRest[idx] == 0.0)
			{
				ratioZero[idx] = 1.0 - wZero/totalZero;
				ratioRest[idx] = 1.0 - wRest/totalRest;
			}
		}

		cout << endl << "\tThresh\t%0>T\t%Spec > T" << endl;
		for (size_t i=0; i<ratioZero.size(); i++)
		{
			cout << i << "\t" << setprecision(2) << i*0.01 << "\t" << setprecision(5)
				<< ratioZero[i] << "\t" << ratioRest[i] << endl;

			if (i>= 10)
				i++;
			if (i>= 20)
				i+=3;
		}
	}
	return threshold;
}



// computes the sqs score using only the regression models
// The computaiton is done on a filtered version of the spectrum (that uses a window of 15 peaks)
double PMCSQS_Scorer::calculateSqsScore(const Config* config, const PeakList& pl, size_t* selectedCharge)
{
	static PeakList sqsPeakList;

	sqsPeakList.copyPeakListLocally(pl);
	sqsPeakList.filterWeakPeaks(config, 0, 15); // use fixed value of 15 peaks per 200 Da, so SQS doesn't change
												// according to the level of weak peak filtering
	const size_t sizeIndex = getSqsSizeIndex(pl.getHeader()->getMOverZ());

	initializeForCurrentSpectrum(config, sqsPeakList);
	calculatePmcValuesForSqs(sqsPeakList, bin_increment);

	MlSample sample;
	sample.label = 0;
	fillSqsSample(sqsPeakList, sample);
	const double sqs = calculateSqsScoreFromSample(sample, sizeIndex);

	if (selectedCharge)
		*selectedCharge=selectChargeForSample(sample, sizeIndex);
	
	return sqs;
}



double PMCSQS_Scorer::calculateSqsScoreFromSample(const MlSample& sample, size_t sizeIndex) const
{
	if (regressionSqsModels_.size()<= sizeIndex)
		error("Not enough sqs regression model for size ",sizeIndex);

	const MlScoreModel* scoreModel = regressionSqsModels_[sizeIndex]->getScoreModel();
	assert(scoreModel);

	double p=1.0 - scoreModel->calcScoreForClass(&sample,0);

	// correct prob according to maxCharge
	const double thresh = sqsRegressionNormalizingConstants_[sizeIndex];
	if (p < thresh)
	{
		assert(thresh>0.0);
		p *= (0.1 / thresh);
	}
	else
		p = 0.1 + 0.9*(p - thresh)/(1.0 - thresh);
	
	return p;
}

int PMCSQS_Scorer::selectChargeForSample(const MlSample& sample, size_t sizeIndex,
										 double* p, int* runnerUp) const
{
	if (p)
		*p=0.0;
	if (runnerUp)
		*runnerUp=0;
	
	vector< vector<double> > probs (maximalChargeWithModels_+1, vector<double>(maximalChargeWithModels_+1,-1.0));
	// compute all pairs
	bool foundModels=false;
	for (int charge1=1; charge1<maximalChargeWithModels_; charge1++)
		for (int charge2=charge1+1; charge2<=maximalChargeWithModels_; charge2++)
		{
			if (! regressionChargeModels_[sizeIndex][charge1][charge2])
				continue;

			const MlScoreModel* scoreModel = regressionChargeModels_[sizeIndex][charge1][charge2]->getScoreModel();
			assert( scoreModel );
			double p=scoreModel->calcScoreForClass(&sample,0);
			probs[charge1][charge2]=p;
			probs[charge2][charge1]=1.0-p;
			foundModels = true;
		}

	if (! foundModels)
		return 0;

	// select charge that has maximal number of wins (if tie, choose the one that has 
	// minimum probs when compared to the others)
	int bestCharge=0;
	int bestNumWins=0;
	double bestSumProbs=0.0;
	for (int charge1=1; charge1<probs.size(); charge1++)
	{
		int numWins=0;
		double sumProbs=0.0;
		for (int charge2=1; charge2<probs.size(); charge2++)
		{
			if (probs[charge1][charge2]<0.0)
				continue;
			sumProbs+=probs[charge1][charge2];
			if (probs[charge1][charge2]>0.5)
				numWins++;
		}
		if (numWins>bestNumWins || numWins == bestNumWins && sumProbs>bestSumProbs)
		{
			bestCharge = charge1;
			bestNumWins = numWins;
			bestSumProbs = sumProbs;
		}
	}

	if (bestCharge <= 0)
		return 0;

	// find the charge that has the highest prob compared to the best charge
	if (p || runnerUp)
	{
		double bestOtherProb=0.0;
		int	   bestOtherCharge=0;
		for (int charge = 1; charge<probs.size(); charge++)
		{
			if (probs[bestCharge][charge]>0.0 && probs[bestCharge][charge]<bestOtherProb)
			{
				bestOtherProb = probs[bestCharge][charge];
				bestOtherCharge = charge;
			}
		}
		if (p)
		{
			// there might be strange cases where the charge selected might actually loose
			// to some other charge (but in the overall case it had a better number of wins)
			// in this case give it prob 0.501
			*p = (bestOtherProb < 0.5 ? bestOtherProb : 0.4999);
		}

		if (runnerUp)
			*runnerUp = bestOtherCharge;
	}
	return bestCharge;
}



bool PMCSQS_Scorer::readSqsModels(const Config* config, const char* fileName)
{
	// TODO resolve SQS tretas
	return false;
	config_ = config;

	string path;
	path = config_->get_resource_dir() + "/" + std::string(fileName);
	
	ifstream in_stream(path.c_str());
	if (! in_stream.good())
	{
		cout << "Warning: couldn't open sqs model for reading: " << path << endl;
		return false;
	}

	char buff[512];
	in_stream.getline(buff,256);
	istringstream iss(buff);

	iss >> maximalChargeWithModels_;
	
	in_stream.getline(buff,256);
	iss.clear();
	iss.str(buff);

	int numSizes=0;
	iss >> numSizes;
	assert(numSizes>0);

	sqsMassThresholds_.resize(numSizes,MAX_FLOAT);
	int i;
	for (i=0; i<numSizes; i++)
		iss >> sqsMassThresholds_[i];

	in_stream.getline(buff,256);
	iss.clear();
	iss.str(buff);

	numSizes=0;
	iss >> numSizes;
	assert(numSizes>0);

	sqsRegressionNormalizingConstants_.resize(numSizes);
	regressionSqsModels_.resize(numSizes);
	regressionChargeModels_.resize(numSizes);
	for (i=0; i<numSizes; i++)
	{
		regressionChargeModels_[i].resize(maximalChargeWithModels_+1);
		for (int j=0; j<=maximalChargeWithModels_; j++)
			regressionChargeModels_[i][j].resize(maximalChargeWithModels_+1, 0);
	}

	for (i=0; i<numSizes; i++)
		iss >> sqsRegressionNormalizingConstants_[i];
	
	while (in_stream.getline(buff,256))
	{
		istringstream iss(buff);
		int sizeIndex=-1,charge1=-1, charge2=-1;
		iss >> sizeIndex >> charge1 >> charge2;
		
		assert( sizeIndex>=0 && charge1 >=0);
		
		if (charge1 == 0)
		{
			assert( sizeIndex < regressionSqsModels_.size());
			if (! regressionSqsModels_[sizeIndex])
				  regressionSqsModels_[sizeIndex] = new MlModel;
			regressionSqsModels_[sizeIndex]->readModel(in_stream, config_->get_resource_dir().c_str());
		}
		else
		{
			if (! regressionChargeModels_[sizeIndex][charge1][charge2])
				regressionChargeModels_[sizeIndex][charge1][charge2] = new MlModel;
			regressionChargeModels_[sizeIndex][charge1][charge2]->readModel(in_stream, config_->get_resource_dir().c_str());
		}
	}

	in_stream.close();
	indInitializedSqs_ = true;
	return true;
}
void PMCSQS_Scorer::writeSqsModels(const char* path) const
{
	ofstream out_stream(path,ios::out);
	if (! out_stream.good())
	{
		cout << "Error: couldn't open pmc model for writing: " << path << endl;
		exit(1);
	}

	out_stream << maximalChargeWithModels_ << endl;
	out_stream << sqsMassThresholds_.size() << setprecision(3) << fixed;
	for (size_t i=0; i<sqsMassThresholds_.size(); i++)
		out_stream << "\t" << sqsMassThresholds_[i];
	out_stream << endl;

	out_stream << sqsRegressionNormalizingConstants_.size();
	for (size_t i=0; i<sqsRegressionNormalizingConstants_.size(); i++)
		out_stream << "\t" << sqsRegressionNormalizingConstants_[i];
	out_stream << endl;
	
	for (size_t i=0; i<regressionSqsModels_.size(); i++)
		if (regressionSqsModels_[i])
		{
			out_stream << i << "\t0" << endl;
			regressionSqsModels_[i]->writeModel(out_stream);
		}

	for (size_t i=0; i<regressionChargeModels_.size(); i++)
		for (size_t j=0; j<regressionChargeModels_[i].size(); j++)
			for (size_t k=0; k<regressionChargeModels_[i][j].size(); k++)
				if (regressionChargeModels_[i][j][k])
				{
					out_stream << i << "\t" << j << "\t" << k << endl;
					regressionChargeModels_[i][j][k]->writeModel(out_stream);
				}

	out_stream.close();
}



/****************************************************************************
Finds the bin which has the optimal values (look for the maximal number of pairs).
Performs search near the peptide's true m/z value to compensate for systematic bias
in the precursor mass.
*****************************************************************************/
size_t PMCSQS_Scorer::findOptimalBinIdx(int true_mz_bin, int charge) const
{
	assert(charge<6);
	const int max_bin_offset = 6-charge; // look in the range +- of this value
	const vector<PMCRankStats>& pmcStatistics = currentSpectrumPmcTables_[charge];
	const int min_bin_idx = (true_mz_bin - max_bin_offset>=0 ? true_mz_bin - max_bin_offset : 0);
	const int max_bin_idx = (true_mz_bin + max_bin_offset>= pmcStatistics.size() ? pmcStatistics.size()-1 :
								true_mz_bin + max_bin_offset);

	if (pmcStatistics[true_mz_bin].numFragmentPairs==0 &&
		pmcStatistics[true_mz_bin].numCharge2FragmentPairs==0)
		return true_mz_bin;
	
	int   optimal_bin_idx=NEG_INF;
	float max_num_pairs=NEG_INF;
	float best_offset=POS_INF;

	if (pmcStatistics[true_mz_bin].numFragmentPairs>=pmcStatistics[true_mz_bin].numCharge2FragmentPairs)
	{
		float max_num_pairs=0;
		int bin_idx;
		for (bin_idx = min_bin_idx; bin_idx<=max_bin_idx; bin_idx++)
			if (pmcStatistics[bin_idx].numFragmentPairs > max_num_pairs)
				max_num_pairs = pmcStatistics[bin_idx].numFragmentPairs;

		// find minimal offset
		for (bin_idx = min_bin_idx; bin_idx<=max_bin_idx; bin_idx++)
			if (pmcStatistics[bin_idx].numFragmentPairs == max_num_pairs &&
				pmcStatistics[bin_idx].meanOffsetPairs < best_offset)
			{
				optimal_bin_idx = bin_idx;
				best_offset = pmcStatistics[bin_idx].meanOffsetPairs;
			}

		return static_cast<size_t>(optimal_bin_idx);
		
	}
	else
	// use the charge 2 fragment pairs
	{
		float max_num_pairs=0; 
		int bin_idx;
		for (bin_idx = min_bin_idx; bin_idx<=max_bin_idx; bin_idx++)
			if (pmcStatistics[bin_idx].numCharge2FragmentPairs > max_num_pairs)
				max_num_pairs = pmcStatistics[bin_idx].numCharge2FragmentPairs;

		// find minimal offset
		for (bin_idx = min_bin_idx; bin_idx<=max_bin_idx; bin_idx++)
			if (pmcStatistics[bin_idx].numCharge2FragmentPairs == max_num_pairs &&
				pmcStatistics[bin_idx].meanOffsetCharge2Pairs < best_offset)
			{
				optimal_bin_idx = bin_idx;
				best_offset = pmcStatistics[bin_idx].meanOffsetCharge2Pairs;
			}

		return static_cast<size_t>(optimal_bin_idx);	
	}


	return MAX_SIZE_T;
}



/*********************************************************************************

**********************************************************************************/
void PMCSQS_Scorer::selectTrainingSampleIndexes(
					int charge,
					const vector<MlSample>& samples,
					const PeakList& pl,
					size_t& correctIndex,
					vector<size_t>& badPmcIndexes) const
{
	const vector<PMCRankStats>& pmcStatistics = currentSpectrumPmcTables_[charge];

	Peptide peptide;
	peptide.parseFromString(config_, pl.getHeader()->getPeptideStr());
	if (peptide.get_num_aas()<=0)
	{
		cout << "Error: supplied training spectrum without peptide!" << endl;
		exit(1);
	}
	peptide.calc_mass(config_);
	const mass_t peptideMass = peptide.get_mass()+MASS_H2O;
	const size_t sizeIndex = calcRankModelSizeIndex(charge, peptideMass);
	const mass_t trueMz = (peptideMass + charge*MASS_PROTON)/static_cast<mass_t>(charge);
	const mass_t observedMz = pl.getHeader()->getMOverZ();

	// check that the training sample has an ok offset
	if (fabs(trueMz-observedMz)>10.0)
	{
		
		cout << "Erorr in m/z offsets (remove this spectrum from training set): " << endl;
		cout << fixed << setprecision(2) << "file m/z: " << observedMz << "\t" << 
			"\"true\" m/z: " << trueMz << "\t peptide: " << peptide.as_string(config_) << endl;
		cout << "spectrum: " << pl.getHeader()->getTitle() << endl;
		
		cout << "Mass Cys = " << config_->get_aa2mass()[Cys] << endl;

		exit(1);
	}

	// find the entry with the correct m/z

	size_t idx=0;
	while (idx<pmcStatistics.size() && pmcStatistics[idx].m_over_z<trueMz)
		idx++;

	if (idx>= pmcStatistics.size())
		idx--;

	if (idx>0 && pmcStatistics[idx].m_over_z-trueMz>trueMz-pmcStatistics[idx-1].m_over_z)
		idx--;

	correctIndex = idx;

	vector<int> idxs;
	idxs.clear();
	badPmcIndexes.clear();


	idxs.push_back(correctIndex+4);
	idxs.push_back(correctIndex+5);
	idxs.push_back(correctIndex+6);
	idxs.push_back(correctIndex+9);
	idxs.push_back(correctIndex+10);
	idxs.push_back(correctIndex+11);
	idxs.push_back(correctIndex+15);
	idxs.push_back(correctIndex+19);
	idxs.push_back(correctIndex+20);
	idxs.push_back(correctIndex-4);
	idxs.push_back(correctIndex-5);
	idxs.push_back(correctIndex-6);
	idxs.push_back(correctIndex-9);
	idxs.push_back(correctIndex-10);
	idxs.push_back(correctIndex-11);
	idxs.push_back(correctIndex-15);
	idxs.push_back(correctIndex-19);
	idxs.push_back(correctIndex-20);

	// select upto 7 random samples (make sure they are not close to the correct one)
	for (size_t i=0; i<7; i++)
	{
		const int idx = static_cast<int>(myRandom()*pmcStatistics.size());
		if (abs(static_cast<int>(correctIndex)-idx)<7)
			continue;

		idxs.push_back(idx);
	}

	sort(idxs.begin(),idxs.end());
	for (size_t i=0; i<idxs.size(); i++)
		if (idxs[i]>=0 && idxs[i]<pmcStatistics.size())
			badPmcIndexes.push_back(static_cast<size_t>(idxs[i]));

}



/*************************************************************************
Tests the performance of precursor mass correction
**************************************************************************/
void PMCSQS_Scorer::test_pmc(Config *config, char *specs_file, int charge, 
							 mass_t minMass, mass_t maxMass)
{
/*	BasicSpecReader bsr;
	static QCPeak peaks[5000];

	FileManager fm;
	FileSet fs;
		
	fm.init_from_file(config,specs_file);
	fs.select_files_in_mz_range(fm,minMass,maxMass,charge);

	const int max_to_read_per_file = 5000;

	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();
	const int num_samples = (all_ssf.size()<max_to_read_per_file ? all_ssf.size() :
									max_to_read_per_file);
	
	vector<mass_t> org_offsets;
	vector<mass_t> corr_offsets;

	vector<size_t> ssf_idxs;
	if (num_samples<all_ssf.size())
	{
		chooseKFromN(num_samples,all_ssf.size(),ssf_idxs);
	}
	else
	{
		int i;
		ssf_idxs.resize(all_ssf.size());
		for (i=0; i<all_ssf.size(); i++)
			ssf_idxs[i]=i;
	}

	vector<SingleSpectrumFile *> ssfs;
	int i;
	for (i=0; i<num_samples; i++)
		ssfs.push_back( all_ssf[ssf_idxs[i]]);
	
	//output_pmc_rank_results(fm,charge,ssfs);*/
}


/***********************************************************************************

Functions for training set.


************************************************************************************/

struct ScanPair {
	ScanPair(int f,int sc, string& se) : file_idx(f), scan(sc), seq(se) {};
	ScanPair(int f,int s) : file_idx(f), scan(s) {};
	ScanPair() : file_idx(-1), scan(-1) {};

	bool operator< (const ScanPair& other) const
	{
		return (file_idx<other.file_idx || 
			    (file_idx == other.file_idx && scan<other.scan));
	}

	bool operator == (const ScanPair& other) const
	{
		return (file_idx == other.file_idx && scan == other.scan);
	}


	int file_idx;
	int scan;
	string seq;
};

void read_idxs_from_file(char *file, vector<ScanPair>& final_pairs, int max_size)
{
	ifstream inp(file,ios::in);
	
	if (! inp.good())
	{	
		cout << "Error opening: " << file << endl;
		exit(1);
	}

	vector<ScanPair> pairs;
	pairs.clear();

	char buff[256];
	while (inp.getline(buff,256))
	{
		istringstream iss(buff);
		int f,s;
		string seq;

		iss >> f >> s >> seq;

		
		if (f>=0 && s>=0)
		{
			if (seq.length()>2)
			{
				pairs.push_back(ScanPair(f,s,seq));
			}
			else
				pairs.push_back(ScanPair(f,s));
		}


	}
	inp.close();

	if (pairs.size() > max_size)
	{
		vector<size_t> idxs;
		chooseKFromN(max_size,pairs.size(),idxs);
		final_pairs.resize(max_size);
		int i;
		for (i=0; i<max_size; i++)
			final_pairs[i]=pairs[idxs[i]];
	}
	else
	{
		final_pairs=pairs;
	}

	sort(final_pairs.begin(),final_pairs.end());
}


void PMCSQS_Scorer::fillRankboostPmcSamplesOld(const PeakList& pl,
								 int charge,
								 vector<RankBoostSample>& samples) const
{
	const int numSamples = currentSpectrumPmcTables_[charge].size();
	const int idx_skip = int((1.0/bin_increment)+0.00001);
	vector<int> idx_offsets;
	int i;

	idx_offsets.clear();
	idx_offsets.push_back(-2*idx_skip);
	idx_offsets.push_back(-1*idx_skip);
	idx_offsets.push_back(idx_skip);
	idx_offsets.push_back(2*idx_skip);

	if (samples.size() != numSamples)
		samples.resize(numSamples);

	for (i=0; i<numSamples; i++)
	{
		const PMCRankStats& stats = currentSpectrumPmcTables_[charge][i];
		RankBoostSample& sam = samples[i];

		const float inten_norm = 1.0/(currentSpectrumTotalIntensity_+1.0);
		int r_idx=0;
		const mass_t mz_offset = (stats.m_over_z - pl.getHeader()->getMOverZ());

		sam.clear();
		sam.add_real_feature(r_idx++,mz_offset);

		if (stats.numFragmentPairs<=2)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else if (stats.numFragmentPairs<4)
		{
			sam.add_real_feature(r_idx+1,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+2,mz_offset);

		r_idx+=3;

		if (stats.numStrongFragmentPairs<3)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+1,mz_offset);

		r_idx+=2;

		if (stats.numCharge2FragmentPairs<=2)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else if (stats.numCharge2FragmentPairs<4)
		{
			sam.add_real_feature(r_idx+1,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+2,mz_offset);

		r_idx+=3;

		if (stats.numStrongCharge2FragmentPairs<3)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+1,mz_offset);

		r_idx+=2;

			
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

		sam.add_real_feature(r_idx++,stats.numFragmentPairs);
		sam.add_real_feature(r_idx++,stats.numStrongFragmentPairs);
		sam.add_real_feature(r_idx++,stats.numCharge2FragmentPairs);
		sam.add_real_feature(r_idx++,stats.numStrongCharge2FragmentPairs);
		sam.add_real_feature(r_idx++,stats.numH2oLossFragmentPairs);
		sam.add_real_feature(r_idx++,stats.numCharge2H2oLossFragmentPairs);

		sam.add_real_feature(r_idx++,stats.intensityInFragmentPairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.intensityInStrongPairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.intensityInCharge2Pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.intensityInStrongCharge2Pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.intensityInH2oLossPairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.intensityInCharge2H2oLossPairs * inten_norm);

		// averages of top k offsets

		float avg=0;
		int j;
		for (j =0; j<7 && j<stats.offsetPairsOderedByIntensity.size(); j++)
		{
			avg += fabs(stats.offsetPairsOderedByIntensity[j]);
			if (j>=2)
				sam.add_real_feature(r_idx+j-2,avg/(float)j);
		}
		r_idx+=5;

		avg=0;
		for (j =0; j<7 && j<stats.charge2offsetPairsOderedByIntensity.size(); j++)
		{
			avg += fabs(stats.charge2offsetPairsOderedByIntensity[j]);
			if (j>=2)
				sam.add_real_feature(r_idx+j-2,avg/(float)j);
		}
		r_idx+=5;


		// offset data
	
		if (stats.meanOffsetPairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.meanOffsetPairs);
			sam.add_real_feature(r_idx++,stats.meanOffsetPairs/(1.0+stats.numFragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.meanOffsetStrongPairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.meanOffsetStrongPairs);
			sam.add_real_feature(r_idx++,stats.meanOffsetStrongPairs/(1.0+stats.numStrongFragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.meanOffsetCharge2Pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.meanOffsetCharge2Pairs);
			sam.add_real_feature(r_idx++,stats.meanOffsetCharge2Pairs/(1.0+stats.numCharge2FragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.meanOffsetCharge2StrongPairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.meanOffsetCharge2StrongPairs);
			sam.add_real_feature(r_idx++,stats.meanOffsetCharge2StrongPairs/(1.0+stats.numStrongCharge2FragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.meanOffsetH2oPairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.meanOffsetH2oPairs);
			sam.add_real_feature(r_idx++,stats.meanOffsetH2oPairs/(1.0+stats.numH2oLossFragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.meanOffsetCharge2H2oPairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.meanOffsetCharge2H2oPairs);
			sam.add_real_feature(r_idx++,stats.meanOffsetCharge2H2oPairs/(1.0+stats.numCharge2H2oLossFragmentPairs));
		}
		else
			r_idx+=2;

		// individual offsets
		for (j=0; j<5 && j<stats.offsetPairsOderedByIntensity.size(); j++)
			sam.add_real_feature(r_idx+j,stats.offsetPairsOderedByIntensity[j]);
		r_idx+=5;

		for (j=0; j<5 && j<stats.charge2offsetPairsOderedByIntensity.size(); j++)
			sam.add_real_feature(r_idx+j,stats.charge2offsetPairsOderedByIntensity[j]);
		r_idx+=5;

	
		// add the +0 +1 +2 strict counts
		sam.add_real_feature(r_idx++,stats.numPairsIso0);
		sam.add_real_feature(r_idx++,stats.intensityInPairsIso0 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.numPairsIso1);
		sam.add_real_feature(r_idx++,stats.intensityInPairsIso1 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.numPairsIso2);
		sam.add_real_feature(r_idx++,stats.intensityInPairsIso2 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.numCharge2PairsIso0);
		sam.add_real_feature(r_idx++,stats.intensityInCharge2PairsIso0 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.numCharge2PairsIso1);
		sam.add_real_feature(r_idx++,stats.intensityInCharge2PairsIso1 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.numCharge2PairsIso2);
		sam.add_real_feature(r_idx++,stats.intensityInCharge2PairsIso2 * inten_norm);

		// add comparative features to -2 -1 +1 +2 Da away
		for (j=0; j<idx_offsets.size(); j++)
		{
			const int other_idx = i + idx_offsets[j];
			if (other_idx<0 || other_idx>= samples.size())
			{
				r_idx+=12;
				continue;
			}

			const PMCRankStats& other = currentSpectrumPmcTables_[charge][other_idx];

			sam.add_real_feature(r_idx++,stats.numFragmentPairs - other.numFragmentPairs);
			sam.add_real_feature(r_idx++,stats.numStrongFragmentPairs - other.numStrongFragmentPairs);
			sam.add_real_feature(r_idx++,stats.numCharge2FragmentPairs - other.numCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,stats.numStrongCharge2FragmentPairs - other.numStrongCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,stats.numH2oLossFragmentPairs - other.numH2oLossFragmentPairs);
			sam.add_real_feature(r_idx++,stats.numCharge2H2oLossFragmentPairs - other.numCharge2H2oLossFragmentPairs);

			sam.add_real_feature(r_idx++,(stats.intensityInFragmentPairs - other.intensityInFragmentPairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.intensityInStrongPairs - other.intensityInStrongPairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.intensityInCharge2Pairs - other.intensityInCharge2Pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.intensityInStrongCharge2Pairs - other.intensityInStrongCharge2Pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.intensityInH2oLossPairs - other.intensityInH2oLossPairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.intensityInCharge2H2oLossPairs - other.intensityInCharge2H2oLossPairs) * inten_norm);
		}

		const int plus_idx = i + idx_skip;
		const int minus_idx = i-idx_skip;

		if (plus_idx<samples.size() && minus_idx>0)
		{
			const PMCRankStats& plus = currentSpectrumPmcTables_[charge][plus_idx];
			const PMCRankStats& minus = currentSpectrumPmcTables_[charge][minus_idx];

			sam.add_real_feature(r_idx++,plus.numFragmentPairs - minus.numFragmentPairs);
			sam.add_real_feature(r_idx++,plus.numStrongFragmentPairs - minus.numStrongFragmentPairs);
			sam.add_real_feature(r_idx++,plus.numCharge2FragmentPairs - minus.numCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,plus.numStrongCharge2FragmentPairs - minus.numStrongCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,plus.numH2oLossFragmentPairs - minus.numH2oLossFragmentPairs);
			sam.add_real_feature(r_idx++,plus.numCharge2H2oLossFragmentPairs - minus.numCharge2H2oLossFragmentPairs);

			sam.add_real_feature(r_idx++,(plus.intensityInFragmentPairs - minus.intensityInFragmentPairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.intensityInStrongPairs - minus.intensityInStrongPairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.intensityInCharge2Pairs - minus.intensityInCharge2Pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.intensityInStrongCharge2Pairs - minus.intensityInStrongCharge2Pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.intensityInH2oLossPairs - minus.intensityInH2oLossPairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.intensityInCharge2H2oLossPairs - minus.intensityInCharge2H2oLossPairs) * inten_norm);
		}
	}
}

/*
void PMCSQS_Scorer::fillRankboostPmcSamplesOld(const PeakList& pl,
								 int charge,
								 vector<RankBoostSample>& samples) const
{
	const int numSamples = currentSpectrumPmcTables_[charge].size();
	const int idx_skip = int((1.0/bin_increment)+0.00001);
	vector<int> idx_offsets;
	int i;

	idx_offsets.clear();
	idx_offsets.push_back(-2*idx_skip);
	idx_offsets.push_back(-1*idx_skip);
	idx_offsets.push_back(idx_skip);
	idx_offsets.push_back(2*idx_skip);

	if (samples.size() != numSamples)
		samples.resize(numSamples);

	for (i=0; i<numSamples; i++)
	{
		const PMCRankStats& stats = currentSpectrumPmcTables_[charge][i];
		RankBoostSample& sam = samples[i];

		const float inten_norm = 1.0/(curr_spec_total_intensity+1.0);
		int r_idx=0;
		const mass_t mz_offset = (stats.m_over_z - pl.getHeader()->getMOverZ());

		sam.clear();
		sam.add_real_feature(r_idx++,mz_offset);

		if (stats.numFragmentPairs<=2)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else if (stats.numFragmentPairs<4)
		{
			sam.add_real_feature(r_idx+1,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+2,mz_offset);

		r_idx+=3;

		if (stats.numStrongFragmentPairs<3)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+1,mz_offset);

		r_idx+=2;

		if (stats.numCharge2FragmentPairs<=2)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else if (stats.numCharge2FragmentPairs<4)
		{
			sam.add_real_feature(r_idx+1,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+2,mz_offset);

		r_idx+=3;

		if (stats.numStrongCharge2FragmentPairs<3)
		{
			sam.add_real_feature(r_idx,mz_offset);
		}
		else
			sam.add_real_feature(r_idx+1,mz_offset);

		r_idx+=2;

			
	//	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=2");
	//	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS <=5");
	//	names.push_back("OFFSET FROM MEASURED M/Z, NUM PAIRS >5");
	//	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS <4");
	//	names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG PAIRS >4");
//
//		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=2");
//		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS <=5");
//		names.push_back("OFFSET FROM MEASURED M/Z, NUM C2 PAIRS >5");
//		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS <4");
//		names.push_back("OFFSET FROM MEASURED M/Z, NUM STRONG C2 PAIRS >4");

		sam.add_real_feature(r_idx++,stats.numFragmentPairs);
		sam.add_real_feature(r_idx++,stats.numStrongFragmentPairs);
		sam.add_real_feature(r_idx++,stats.numCharge2FragmentPairs);
		sam.add_real_feature(r_idx++,stats.numStrongCharge2FragmentPairs);
		sam.add_real_feature(r_idx++,stats.numH2oLossFragmentPairs);
		sam.add_real_feature(r_idx++,stats.numCharge2H2oLossFragmentPairs);

		sam.add_real_feature(r_idx++,stats.intensityInFragmentPairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.intensityInStrongPairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.intensityInCharge2Pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.intensityInStrongCharge2Pairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.intensityInH2oLossPairs * inten_norm);
		sam.add_real_feature(r_idx++,stats.intensityInCharge2H2oLossPairs * inten_norm);

		// averages of top k offsets

		float avg=0;
		int j;
		for (j =0; j<7 && j<stats.offsetPairsOderedByIntensity.size(); j++)
		{
			avg += fabs(stats.offsetPairsOderedByIntensity[j]);
			if (j>=2)
				sam.add_real_feature(r_idx+j-2,avg/(float)j);
		}
		r_idx+=5;

		avg=0;
		for (j =0; j<7 && j<stats.charge2offsetPairsOderedByIntensity.size(); j++)
		{
			avg += fabs(stats.charge2offsetPairsOderedByIntensity[j]);
			if (j>=2)
				sam.add_real_feature(r_idx+j-2,avg/(float)j);
		}
		r_idx+=5;


		// offset data
	
		if (stats.meanOffsetPairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.meanOffsetPairs);
			sam.add_real_feature(r_idx++,stats.meanOffsetPairs/(1.0+stats.numFragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.meanOffsetStrongPairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.meanOffsetStrongPairs);
			sam.add_real_feature(r_idx++,stats.meanOffsetStrongPairs/(1.0+stats.numStrongFragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.meanOffsetCharge2Pairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.meanOffsetCharge2Pairs);
			sam.add_real_feature(r_idx++,stats.meanOffsetCharge2Pairs/(1.0+stats.numCharge2FragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.meanOffsetCharge2StrongPairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.meanOffsetCharge2StrongPairs);
			sam.add_real_feature(r_idx++,stats.meanOffsetCharge2StrongPairs/(1.0+stats.numStrongCharge2FragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.meanOffsetH2oPairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.meanOffsetH2oPairs);
			sam.add_real_feature(r_idx++,stats.meanOffsetH2oPairs/(1.0+stats.numH2oLossFragmentPairs));
		}
		else
			r_idx+=2;

		if (stats.meanOffsetCharge2H2oPairs<POS_INF)
		{
			sam.add_real_feature(r_idx++,stats.meanOffsetCharge2H2oPairs);
			sam.add_real_feature(r_idx++,stats.meanOffsetCharge2H2oPairs/(1.0+stats.numCharge2H2oLossFragmentPairs));
		}
		else
			r_idx+=2;

		// individual offsets
		for (j=0; j<5 && j<stats.offsetPairsOderedByIntensity.size(); j++)
			sam.add_real_feature(r_idx+j,stats.offsetPairsOderedByIntensity[j]);
		r_idx+=5;

		for (j=0; j<5 && j<stats.charge2offsetPairsOderedByIntensity.size(); j++)
			sam.add_real_feature(r_idx+j,stats.charge2offsetPairsOderedByIntensity[j]);
		r_idx+=5;

	
		// add the +0 +1 +2 strict counts
		sam.add_real_feature(r_idx++,stats.numPairsIso0);
		sam.add_real_feature(r_idx++,stats.intensityInPairsIso0 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.numPairsIso1);
		sam.add_real_feature(r_idx++,stats.intensityInPairsIso1 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.numPairsIso2);
		sam.add_real_feature(r_idx++,stats.intensityInPairsIso2 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.numCharge2PairsIso0);
		sam.add_real_feature(r_idx++,stats.intensityInCharge2PairsIso0 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.numCharge2PairsIso1);
		sam.add_real_feature(r_idx++,stats.intensityInCharge2PairsIso1 * inten_norm);
	
		sam.add_real_feature(r_idx++,stats.numCharge2PairsIso2);
		sam.add_real_feature(r_idx++,stats.intensityInCharge2PairsIso2 * inten_norm);

		// add comparative features to -2 -1 +1 +2 Da away
		for (j=0; j<idx_offsets.size(); j++)
		{
			const int other_idx = i + idx_offsets[j];
			if (other_idx<0 || other_idx>= samples.size())
			{
				r_idx+=12;
				continue;
			}

			const PMCRankStats& other = currentSpectrumPmcTables_[charge][other_idx];

			sam.add_real_feature(r_idx++,stats.numFragmentPairs - other.numFragmentPairs);
			sam.add_real_feature(r_idx++,stats.numStrongFragmentPairs - other.numStrongFragmentPairs);
			sam.add_real_feature(r_idx++,stats.numCharge2FragmentPairs - other.numCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,stats.numStrongCharge2FragmentPairs - other.numStrongCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,stats.numH2oLossFragmentPairs - other.numH2oLossFragmentPairs);
			sam.add_real_feature(r_idx++,stats.numCharge2H2oLossFragmentPairs - other.numCharge2H2oLossFragmentPairs);

			sam.add_real_feature(r_idx++,(stats.intensityInFragmentPairs - other.intensityInFragmentPairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.intensityInStrongPairs - other.intensityInStrongPairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.intensityInCharge2Pairs - other.intensityInCharge2Pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.intensityInStrongCharge2Pairs - other.intensityInStrongCharge2Pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.intensityInH2oLossPairs - other.intensityInH2oLossPairs) * inten_norm);
			sam.add_real_feature(r_idx++,(stats.intensityInCharge2H2oLossPairs - other.intensityInCharge2H2oLossPairs) * inten_norm);
		}

		const int plus_idx = i + idx_skip;
		const int minus_idx = i-idx_skip;

		if (plus_idx<samples.size() && minus_idx>0)
		{
			const PMCRankStats& plus = currentSpectrumPmcTables_[charge][plus_idx];
			const PMCRankStats& minus = currentSpectrumPmcTables_[charge][minus_idx];

			sam.add_real_feature(r_idx++,plus.numFragmentPairs - minus.numFragmentPairs);
			sam.add_real_feature(r_idx++,plus.numStrongFragmentPairs - minus.numStrongFragmentPairs);
			sam.add_real_feature(r_idx++,plus.numCharge2FragmentPairs - minus.numCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,plus.numStrongCharge2FragmentPairs - minus.numStrongCharge2FragmentPairs);
			sam.add_real_feature(r_idx++,plus.numH2oLossFragmentPairs - minus.numH2oLossFragmentPairs);
			sam.add_real_feature(r_idx++,plus.numCharge2H2oLossFragmentPairs - minus.numCharge2H2oLossFragmentPairs);

			sam.add_real_feature(r_idx++,(plus.intensityInFragmentPairs - minus.intensityInFragmentPairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.intensityInStrongPairs - minus.intensityInStrongPairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.intensityInCharge2Pairs - minus.intensityInCharge2Pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.intensityInStrongCharge2Pairs - minus.intensityInStrongCharge2Pairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.intensityInH2oLossPairs - minus.intensityInH2oLossPairs) * inten_norm);
			sam.add_real_feature(r_idx++,(plus.intensityInCharge2H2oLossPairs - minus.intensityInCharge2H2oLossPairs) * inten_norm);
		}
	}
}

*/


void create_training_files(Config *config)
{
/*	char mzxml_list[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\HEK293_mzxml_list.txt"};
	char idxs_neg_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\H40ul_neg_samples.txt"};
//	char idxs1_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\H40ul_pos_samples.1.txt"};
//	char idxs2_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\H40ul_pos_samples.2.txt"};
//	char idxs2_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\Len10_pos_samples.2.txt"};
	char idxs1_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\sqs_train_pos_samples.1.txt"};
	char idxs2_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\sqs_train_pos_samples.2.txt"};
	char idxs3_file[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\H40ul_pos_samples.3.txt"};

	char out_base[]={"C:\\Work\\msms5\\PepNovoHQ\\pmcsqs\\sqs_train"};
	string out_neg (out_base); 
	string out1=out_neg;
	string out2=out_neg;
	string out3=out_neg;

	out_neg += "_neg.mgf";
	out1 += "_1.mgf";
	out2 += "_2.mgf";
	out3 += "_3.mgf";

	ofstream stream_neg (out_neg.c_str(),ios::out);
	ofstream stream1(out1.c_str(),ios::out);
	ofstream stream2(out2.c_str(),ios::out);
	ofstream stream3(out3.c_str(),ios::out);

	vector<ScanPair> neg_pairs, pairs1,pairs2,pairs3;


	read_idxs_from_file(idxs_neg_file,neg_pairs,12000);
	read_idxs_from_file(idxs1_file,pairs1,12000);
	read_idxs_from_file(idxs2_file,pairs2,12000);
	read_idxs_from_file(idxs3_file,pairs3,8000);

	cout << "Read " << neg_pairs.size() << " neg idxs\n";
	cout << "Read " << pairs1.size() << " pos1 idxs\n";
	cout << "Read " << pairs2.size() << " pos2 idxs\n";
	cout << "Read " << pairs3.size() << " pos3 idxs\n";

	vector<bool> file_inds;
	file_inds.resize(10000,false);
	int i;

	for (i=0; i<neg_pairs.size(); i++)
		file_inds[neg_pairs[i].file_idx]=true;

	for (i=0; i<pairs1.size(); i++)
		file_inds[pairs1[i].file_idx]=true;

	for (i=0; i<pairs2.size(); i++)
		file_inds[pairs2[i].file_idx]=true;

	for (i=0; i<pairs3.size(); i++)
		file_inds[pairs3[i].file_idx]=true;

	
	FileManager fm;
	FileSet fs;

	fm.init_from_list_file(config,mzxml_list,file_inds);
	fs.select_all_files(fm);
	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();

	

	// read spectra
	BasicSpecReader bsr;
	QCPeak peaks[5000];

	int num_out_neg=0, num_out1=0, num_out2=0, num_out3=0;
	int neg_idx=0,c1_idx=0,c2_idx=0,c3_idx=0;


	for (i=0; i<all_ssf.size(); i++)
	{
		MZXML_single *ssf = (MZXML_single *)all_ssf[i];
		ScanPair ssf_pair(ssf->file_idx,ssf->scan_number);
		string seq="";

		int out_dest=-1;

		while (neg_idx<neg_pairs.size() && neg_pairs[neg_idx]<ssf_pair)
			neg_idx++;

		if (neg_idx<neg_pairs.size() && neg_pairs[neg_idx]==ssf_pair)
			out_dest=0;


		while (c1_idx<pairs1.size() && pairs1[c1_idx]<ssf_pair)
			c1_idx++;
		if (c1_idx<pairs1.size() && pairs1[c1_idx]==ssf_pair)
		{
			seq = pairs1[c1_idx].seq;
			out_dest=1;
		}


		while (c2_idx<pairs2.size() && pairs2[c2_idx]<ssf_pair)
			c2_idx++;
		if (c2_idx<pairs2.size() && pairs2[c2_idx]==ssf_pair)
		{
			seq = pairs2[c2_idx].seq;
			out_dest=2;
		}


		while (c3_idx<pairs3.size() && pairs3[c3_idx]<ssf_pair)
			c3_idx++;
		if (c3_idx<pairs3.size() && pairs3[c3_idx]==ssf_pair)
		{
			seq = pairs3[c3_idx].seq;
			out_dest=3;
		}

		if (out_dest<0)
			continue;

		BasicSpectrum bs;
		bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks);
		bs.peaks = peaks;
		bs.ssf = ssf;

	//	if (out_dest>0)
	//		bs.ssf->peptide.parseFromString(config,seq);
	
		char name_buff[64];
		if (out_dest==0)
		{
			sprintf(name_buff,"train_neg_%d_%d_%d",num_out_neg,ssf->file_idx,ssf->scan_number);
			bs.ssf->single_name = string(name_buff);
			bs.output_to_mgf(stream_neg,config);
			num_out_neg++;
			continue;
		}

		if (out_dest==1)
		{
			sprintf(name_buff,"train_pos1_%d_%d_%d",num_out1,ssf->file_idx,ssf->scan_number);
			bs.ssf->single_name = string(name_buff);
			bs.output_to_mgf(stream1,config,seq.c_str());
			num_out1++;
			continue;
		}

		if (out_dest==2)
		{
			sprintf(name_buff,"train_pos2_%d_%d_%d",num_out2,ssf->file_idx,ssf->scan_number);
			bs.ssf->single_name = string(name_buff);
			bs.output_to_mgf(stream2,config,seq.c_str());
			num_out2++;
			continue;
		}

		if (out_dest==3)
		{
			sprintf(name_buff,"train_pos3_%d_%d_%d",num_out3,ssf->file_idx,ssf->scan_number);
			bs.ssf->single_name = string(name_buff);
			bs.output_to_mgf(stream3,config,seq.c_str());
			num_out3++;
			continue;
		}
	}

	cout << "Wrote: " << endl;
	cout << "Neg " << num_out_neg << " / " << neg_pairs.size() << endl;
	cout << "Pos1 " << num_out1 << " / " << pairs1.size() << endl;
	cout << "Pos2 " << num_out2 << " / " << pairs2.size() << endl;
	cout << "Pos3 " << num_out3 << " / " << pairs3.size() << endl;

	stream_neg.close();
	stream1.close();
	stream2.close();
	stream3.close();*/
	
}










